//! Append-only (write-once) file persistence for an [`AuditTrail`].
//!
//! A [`WormLog`] is a journal file with WORM semantics: it is created once,
//! appended to line by line (each append `fsync`ed before returning), and
//! never rewritten. Every record carries the running hash-chain digest of
//! the entry it stores, so reopening the journal verifies the whole chain
//! and detects retroactive edits, spliced lines from another run, and torn
//! final writes (a crash mid-append leaves a line without its newline).
//!
//! The journal covers the **entries** of a trail; electronic signatures and
//! policy settings are trail-side state and are re-applied after a reload
//! ([`WormLog::into_trail`] returns a fresh, verified trail).

use std::fs::{File, OpenOptions};
use std::io::Write;
use std::path::{Path, PathBuf};

use tpt_med_audit::{hex, sha256};
use tpt_med_core::{AuditAction, AuditEvent};

use crate::{AuditEntry, AuditTrail, UtcStamp};

/// Journal header magic, versioned so a format change is detectable.
const HEADER_MAGIC: &str = "TPT-WORM-1";

/// WORM journal failure modes.
#[derive(Debug)]
pub enum WormError {
    /// Underlying file system error.
    Io(std::io::Error),
    /// The file does not end with a newline: the last append was torn by a
    /// crash. The journal is refused rather than silently truncated.
    TornTail,
    /// The header line is not [`HEADER_MAGIC`] (wrong file or format).
    BadHeader,
    /// A record line could not be parsed.
    MalformedLine {
        /// 1-based line number (header is line 1).
        line: usize,
    },
    /// Sequence numbers are not contiguous.
    SequenceDiscontinuity {
        /// The sequence the journal required.
        expected: u64,
        /// The sequence the record carried.
        found: u64,
    },
    /// The stored chain digest does not match the recomputation — the line
    /// (or an earlier one) was edited after the fact.
    ChainBroken {
        /// 1-based line number of the first broken link.
        line: usize,
    },
    /// A record was spliced in from a different run.
    RunIdMismatch {
        /// The journal's run id.
        expected: String,
        /// The spliced record's own run binding, when detectable.
        found: String,
    },
}

impl core::fmt::Display for WormError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            WormError::Io(e) => write!(f, "io error: {e}"),
            WormError::TornTail => {
                write!(f, "journal ends mid-line (torn final append)")
            }
            WormError::BadHeader => write!(f, "not a {HEADER_MAGIC} journal"),
            WormError::MalformedLine { line } => write!(f, "malformed record at line {line}"),
            WormError::SequenceDiscontinuity { expected, found } => {
                write!(f, "sequence {found} where {expected} required")
            }
            WormError::ChainBroken { line } => {
                write!(f, "chain digest mismatch at line {line}")
            }
            WormError::RunIdMismatch { expected, found } => {
                write!(
                    f,
                    "record from run '{found}' in journal for run '{expected}'"
                )
            }
        }
    }
}

impl std::error::Error for WormError {}

impl From<std::io::Error> for WormError {
    fn from(e: std::io::Error) -> Self {
        WormError::Io(e)
    }
}

/// Journal escaping: the record format is pipe-separated, so `|` and `\`
/// (plus line breaks) are escaped in the free-text-ish fields. The digest,
/// sequence, timestamp and action fields are emitted from their own
/// constrained alphabets.
fn esc(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '|' => out.push_str("\\p"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            c => out.push(c),
        }
    }
    out
}

/// Inverse of [`esc`].
fn unesc(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    let mut chars = s.chars();
    while let Some(c) = chars.next() {
        if c == '\\' {
            match chars.next() {
                Some('\\') => out.push('\\'),
                Some('p') => out.push('|'),
                Some('n') => out.push('\n'),
                Some('r') => out.push('\r'),
                Some(other) => {
                    out.push('\\');
                    out.push(other);
                }
                None => out.push('\\'),
            }
        } else {
            out.push(c);
        }
    }
    out
}

/// An append-only, write-once journal of audit entries.
///
/// Record format (one line each, after a `HEADER_MAGIC|run_id` header):
/// `digest|sequence|epoch_seconds|nanos|actor|object_class|object_token|action|reason`,
/// where `digest` is the running hash-chain link `H(prev ‖ canonical)` over
/// the entry's [`AuditEntry::canonical`] form — the same chain the in-memory
/// trail verifies, extended one link per journal append.
pub struct WormLog {
    path: PathBuf,
    file: File,
    run_id: String,
    entries: Vec<AuditEntry>,
    prev_digest: [u8; 32],
}

impl core::fmt::Debug for WormLog {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("WormLog")
            .field("path", &self.path)
            .field("run_id", &self.run_id)
            .field("entries", &self.entries.len())
            .finish_non_exhaustive()
    }
}

impl WormLog {
    /// Creates a new journal file. Fails if the path already exists — a
    /// journal is write-once by construction, never reopened for a second
    /// lifetime.
    pub fn create(path: impl AsRef<Path>, run_id: impl Into<String>) -> Result<Self, WormError> {
        let run_id = run_id.into();
        let mut file = OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(path.as_ref())?;
        writeln!(file, "{HEADER_MAGIC}|{run_id}")?;
        file.sync_all()?;
        let prev_digest = sha256(run_id.as_bytes());
        Ok(Self {
            path: path.as_ref().to_path_buf(),
            file,
            run_id,
            entries: Vec::new(),
            prev_digest,
        })
    }

    /// Opens an existing journal, verifying the header, per-line chain
    /// digests, and sequence continuity. A torn final line (no trailing
    /// newline) is refused — the caller decides the forensic response; this
    /// layer does not silently truncate a regulated record.
    pub fn open(path: impl AsRef<Path>) -> Result<Self, WormError> {
        let path = path.as_ref().to_path_buf();
        let raw = std::fs::read(&path)?;
        let text = String::from_utf8(raw).map_err(|_| WormError::BadHeader)?;
        let mut lines = text.split('\n').enumerate().peekable();
        let Some((_, header)) = lines.next() else {
            return Err(WormError::BadHeader);
        };
        let mut header_parts = header.splitn(2, '|');
        if header_parts.next() != Some(HEADER_MAGIC) {
            return Err(WormError::BadHeader);
        }
        let run_id = header_parts.next().ok_or(WormError::BadHeader)?.to_string();

        let mut entries = Vec::new();
        let mut prev_digest = sha256(run_id.as_bytes());
        // A well-formed file ends with '\n', so the final peekable element
        // is the empty string after the split.
        let torn = !text.ends_with('\n');
        while let Some((i, line)) = lines.next() {
            let is_last = lines.peek().is_none();
            if is_last && line.is_empty() && !torn {
                break;
            }
            if is_last && torn {
                return Err(WormError::TornTail);
            }
            let lineno = i + 1;
            let fields: Vec<&str> = line.splitn(9, '|').collect();
            let [digest, seq, secs, nanos, actor, class, token, action, reason] = fields[..] else {
                return Err(WormError::MalformedLine { line: lineno });
            };
            let sequence: u64 = seq
                .parse()
                .map_err(|_| WormError::MalformedLine { line: lineno })?;
            if sequence != entries.len() as u64 {
                return Err(WormError::SequenceDiscontinuity {
                    expected: entries.len() as u64,
                    found: sequence,
                });
            }
            let epoch_seconds: i64 = secs
                .parse()
                .map_err(|_| WormError::MalformedLine { line: lineno })?;
            let nanos: u32 = nanos
                .parse()
                .map_err(|_| WormError::MalformedLine { line: lineno })?;
            let action = parse_action(action).ok_or(WormError::MalformedLine { line: lineno })?;
            let entry = AuditEntry {
                sequence,
                timestamp: UtcStamp {
                    epoch_seconds,
                    nanos,
                },
                actor: unesc(actor),
                object_class: unesc(class),
                object_token: unesc(token),
                action,
                reason: unesc(reason),
            };
            // Chain verification: the recomputed link over the canonical
            // form must equal the stored digest, and the record must bind
            // to this journal's run (a spliced line breaks the chain; the
            // explicit header makes the diagnosis clearer).
            if digest.len() != 64 || !digest.bytes().all(|b| b.is_ascii_hexdigit()) {
                return Err(WormError::MalformedLine { line: lineno });
            }
            let canonical = entry.canonical();
            let mut buf = prev_digest.to_vec();
            buf.extend_from_slice(canonical.as_bytes());
            let expected = hex(&sha256(&buf));
            if expected != digest {
                return Err(WormError::ChainBroken { line: lineno });
            }
            prev_digest = sha256(&buf);
            entries.push(entry);
        }

        let file = OpenOptions::new().append(true).open(&path)?;
        Ok(Self {
            path,
            file,
            run_id,
            entries,
            prev_digest,
        })
    }

    /// Appends one entry (sequence assigned from the journal length), with
    /// the record flushed and `fsync`ed before returning.
    pub fn append(&mut self, event: AuditEvent) -> Result<(), WormError> {
        let entry = AuditEntry {
            sequence: self.entries.len() as u64,
            timestamp: UtcStamp::now(),
            actor: event.actor_token,
            object_class: event.object_class,
            object_token: event.object_token,
            action: event.action,
            reason: event.reason,
        };
        let canonical = entry.canonical();
        let mut buf = self.prev_digest.to_vec();
        buf.extend_from_slice(canonical.as_bytes());
        let digest = hex(&sha256(&buf));
        let line = format!(
            "{digest}|{}|{}|{}|{}|{}|{}|{}|{}\n",
            entry.sequence,
            entry.timestamp.epoch_seconds,
            entry.timestamp.nanos,
            esc(&entry.actor),
            esc(&entry.object_class),
            esc(&entry.object_token),
            entry.action,
            esc(&entry.reason),
        );
        self.file.write_all(line.as_bytes())?;
        self.file.flush()?;
        self.file.sync_all()?;
        self.prev_digest = sha256(&buf);
        self.entries.push(entry);
        Ok(())
    }

    /// The journal's file path.
    pub fn path(&self) -> &Path {
        &self.path
    }

    /// The run id the journal is bound to.
    pub fn run_id(&self) -> &str {
        &self.run_id
    }

    /// The entries loaded so far.
    pub fn entries(&self) -> &[AuditEntry] {
        &self.entries
    }

    /// Entry count.
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// True when no entries have been appended.
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Rebuilds a live [`AuditTrail`] from the journal's entries. The
    /// returned trail re-verifies its own chain; signatures and policies
    /// are trail-side state and are re-applied by the caller.
    pub fn into_trail(self) -> AuditTrail {
        let mut trail = AuditTrail::from_entries(self.run_id, self.entries);
        trail.digest_index = trail.compute_chain();
        trail
    }
}

/// Parses the action token back to an [`AuditAction`].
fn parse_action(s: &str) -> Option<AuditAction> {
    match s {
        "create" => Some(AuditAction::Create),
        "modify" => Some(AuditAction::Modify),
        "delete" => Some(AuditAction::Delete),
        "approve" => Some(AuditAction::Approve),
        "reject" => Some(AuditAction::Reject),
        "export" => Some(AuditAction::Export),
        "simulate" => Some(AuditAction::Simulate),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn event(action: AuditAction, token: &str) -> AuditEvent {
        AuditEvent::new(
            "operator:op1",
            "simulation",
            token,
            action,
            "screening rerun",
        )
    }

    fn tmp_path(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        p.push(format!("tpt-worm-test-{}-{}.log", std::process::id(), name));
        let _ = std::fs::remove_file(&p);
        p
    }

    #[test]
    fn journal_round_trips_into_a_verifying_trail() {
        let path = tmp_path("roundtrip");
        let mut log = WormLog::create(&path, "run-rt").expect("create");
        log.append(event(AuditAction::Create, "s1")).unwrap();
        log.append(event(AuditAction::Simulate, "s1")).unwrap();
        log.append(event(AuditAction::Modify, "s1")).unwrap();
        assert_eq!(log.len(), 3);
        let trail = log.into_trail();
        assert_eq!(trail.run_id, "run-rt");
        assert!(trail.verify_integrity(), "rebuilt trail must verify");
        assert_eq!(trail.len(), 3);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn journal_survives_a_process_restart() {
        let path = tmp_path("restart");
        {
            let mut log = WormLog::create(&path, "run-restart").expect("create");
            log.append(event(AuditAction::Create, "s1")).unwrap();
            log.append(event(AuditAction::Export, "s1")).unwrap();
        }
        // Process restart: reopen from disk only.
        let mut log = WormLog::open(&path).expect("reopen");
        assert_eq!(log.len(), 2);
        assert_eq!(log.run_id(), "run-restart");
        // Sequence continues where the file left off.
        log.append(event(AuditAction::Approve, "s1")).unwrap();
        drop(log);
        let log = WormLog::open(&path).expect("reopen 2");
        assert_eq!(log.len(), 3);
        assert_eq!(log.entries()[2].sequence, 2);
        assert_eq!(log.entries()[2].action, AuditAction::Approve);
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn retroactive_edits_break_the_chain_on_reopen() {
        let path = tmp_path("tamper");
        {
            let mut log = WormLog::create(&path, "run-tamper").expect("create");
            log.append(event(AuditAction::Create, "s1")).unwrap();
            log.append(event(AuditAction::Modify, "s1")).unwrap();
            log.append(event(AuditAction::Export, "s1")).unwrap();
        }
        // Rewrite one character inside the last record's reason field.
        let raw = std::fs::read(&path).unwrap();
        let text = String::from_utf8(raw).unwrap();
        let mut lines: Vec<String> = text.split('\n').map(String::from).collect();
        let last_record = lines.len() - 2; // trailing empty element after the final newline
        assert!(lines[last_record].contains("screening rerun"));
        lines[last_record] = lines[last_record].replace("screening rerun", "screening Xerun");
        let tampered = lines.join("\n");
        std::fs::write(&path, tampered).unwrap();
        match WormLog::open(&path) {
            Err(WormError::ChainBroken { line }) => assert_eq!(line, 4, "first broken link"),
            other => panic!("expected ChainBroken, got {other:?}"),
        }
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn a_torn_final_line_is_refused() {
        let path = tmp_path("torn");
        {
            let mut log = WormLog::create(&path, "run-torn").expect("create");
            log.append(event(AuditAction::Create, "s1")).unwrap();
        }
        // Simulate a crash mid-append: chop the trailing newline.
        let mut raw = std::fs::read(&path).unwrap();
        assert_eq!(raw.pop(), Some(b'\n'));
        std::fs::write(&path, &raw).unwrap();
        assert!(matches!(WormLog::open(&path), Err(WormError::TornTail)));
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn create_refuses_an_existing_journal() {
        let path = tmp_path("exists");
        let _log = WormLog::create(&path, "run-a").expect("create");
        let second = WormLog::create(&path, "run-b");
        assert!(
            matches!(second, Err(WormError::Io(ref e)) if e.kind() == std::io::ErrorKind::AlreadyExists),
            "{second:?}"
        );
        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn spliced_records_from_another_run_are_caught() {
        let path_a = tmp_path("splice-a");
        let path_b = tmp_path("splice-b");
        {
            let mut a = WormLog::create(&path_a, "run-A").expect("create a");
            a.append(event(AuditAction::Create, "s1")).unwrap();
            let mut b = WormLog::create(&path_b, "run-B").expect("create b");
            b.append(event(AuditAction::Create, "s1")).unwrap();
        }
        // Splice run A's record into run B's journal, renumbered to the
        // next sequence so it passes the sequence check and reaches the
        // chain verification.
        let a = std::fs::read_to_string(&path_a).unwrap();
        let record = a.lines().nth(1).unwrap();
        let mut fields: Vec<&str> = record.splitn(9, '|').collect();
        assert_eq!(fields[1], "0");
        fields[1] = "1";
        let spliced = fields.join("|");
        let mut b = std::fs::read_to_string(&path_b).unwrap();
        b.push_str(&spliced);
        b.push('\n');
        std::fs::write(&path_b, b).unwrap();
        // The chain is seeded by the run id, so the spliced link breaks —
        // reported against the spliced line.
        match WormLog::open(&path_b) {
            Err(WormError::ChainBroken { line }) => assert_eq!(line, 3),
            other => panic!("expected ChainBroken, got {other:?}"),
        }
        let _ = std::fs::remove_file(&path_a);
        let _ = std::fs::remove_file(&path_b);
    }

    #[test]
    fn free_text_fields_round_trip_through_the_journal_escaping() {
        let path = tmp_path("escaping");
        let mut log = WormLog::create(&path, "run-esc").expect("create");
        log.append(AuditEvent::new(
            "operator|op\\1",
            "class|two",
            "tok",
            AuditAction::Modify,
            "reason with | pipes \\ and\nnewlines",
        ))
        .unwrap();
        drop(log);
        let log = WormLog::open(&path).expect("reopen");
        let e = &log.entries()[0];
        assert_eq!(e.actor, "operator|op\\1");
        assert_eq!(e.object_class, "class|two");
        assert_eq!(e.reason, "reason with | pipes \\ and\nnewlines");
        let trail = log.into_trail();
        assert!(trail.verify_integrity());
        let _ = std::fs::remove_file(&path);
    }
}
