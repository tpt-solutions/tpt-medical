# RFC 0010: JPEG 2000 Part 2 Multi-Component — Spike Findings

- **Status:** Draft (spike complete; no implementation proposed)
- **Started:** 2026-09-27
- **Crates:** `tpt-med-dicom`

## Summary

A read-only investigation into the one open `tpt-med-dicom` ingestion item:
decoding a genuine JPEG 2000 Part 2 extended multi-component transform. No code
was written. The conclusion is that **the gap is real, no pure-Rust route
exists today, and the cost is dominated by something other than the transform
itself** — namely access to the ISO/IEC 15444-2 specification. The
recommendation is to keep the item deferred, with the reasoning recorded so it
is not re-investigated from scratch, and with the one input that would change
the answer (purchasing the standard) named explicitly.

## What "still rejected" actually means

The claim in `tpt-med-dicom`'s docs is that a real Part 2 transform is refused
rather than mis-decoded. Verified directly in the pinned source,
`pdfluent-jpeg2000` 0.4.0, in both of its marker loops:

- `src/j2c/codestream.rs:111` — the main-header loop matches an explicit list
  (`SIZ`, `COD`, `COC`, `QCD`, `QCC`, `RGN`, `TLM`, `COM`, `PPM`, `CRG`, and the
  no-parameter `0x30..=0x3F` range) and `bail!(MarkerError::Unsupported)` on
  anything else.
- `src/j2c/tile.rs:279` — the same catch-all in the tile-part loop.

Part 2 signals its component transform with the **MCC** (multiple component
collection), **MCT** (multiple component transform) and **MCOD** markers. None
are in either list, so such a codestream is rejected with a named error. The
existing documentation is accurate as written.

Worth being precise about what *is* supported, because the crate is easy to
misread: `src/j2c/mct.rs` exists and implements a multi-component transform —
but only **Part 1's** fixed 3-component RCT/ICT (Annex G.2/G.3), the reversible
and irreversible colour transforms. That is a different feature from Part 2's
arbitrary N-component transform, and its presence is not evidence of Part 2

## The dependency landscape

| Option | Part 2 support | Cost |
|---|---|---|
| `pdfluent-jpeg2000` 0.4.0 (current) | No. Upstream documents "some color spaces from the extensions (15444-2)" — i.e. the JP2 *container* colour handling, not codestream MCT markers. | Would need an upstream change. |
| `oxideav-jpeg2000` | No, and explicit about it: its README states the Part 2 "scaling based" extended-RGN is outside its Part 1 scope and "surfaces a clean error rather than mis-decoding". A careful, byte-exact Part 1 decoder that declines Part 2 rather than guessing. | Not an alternative. |
| OpenJPEG (`libopenjp2`) | **Yes** — full Part 2, including arbitrary MCT. The reference implementation. | See the constraints below. |
| Write it here | Possible in principle. | See the constraints below. |

The two remaining options are both blocked by workspace constraints rather than
by difficulty:

- **OpenJPEG is a C library.** The workspace sets `unsafe_code = "forbid"`
  globally, and `tpt-med-wasm` — a workspace member that depends on
  `tpt-med-dicom` — targets `wasm32`. A `libopenjp2` dependency would need a C
  toolchain at build time and a port to WASM, which is a substantial project in
  its own right. Adopting it would be a change to the workspace's dependency
  posture, not a crate-local decision. It is also worth noting that switching
  codecs to chase Part 2 would discard the existing, verified signed-component
  fix in `jpeg2000.rs` (the DC level-shift workaround) and re-open that
  question against a different implementation.
- **Writing the transform here is blocked on the specification.** The transform
  mathematics is the easy part — an N×M matrix multiply applied after inverse
  quantisation and the inverse wavelet. The hard part is a
  specification-faithful implementation of the MCC/MCT/MCOD marker syntax and
  its edge cases, and ISO/IEC 15444-2 is a paywalled standard. This repository
  contains **no** copy of the JPEG 2000 standards (`docs/` holds only `api`,
  `book`, `rfc`). For comparison, `oxideav-jpeg2000` documents that it was
  written clean-room from local copies of T.800 / 15444-1 — that was a
  precondition for its quality, and we do not hold the equivalent for Part 2.

support.


## Recommendation

**Keep the item deferred, and record why.** Three reasons, in order of weight:

1. **The input that gates it is a purchase, not engineering.** Obtaining
   ISO/IEC 15444-2 is a procurement decision with a budget and a review
   process. Until that is answered, no amount of engineering time moves this
   item, and a spike repeated later will reach the same conclusion.
2. **The clinical prevalence is low.** `.92`/`.93` are *transfer syntax* UIDs;
   a Part 2 file that uses no extended transform is simply a Part 1 codestream
   wearing a different UID, and those decode today. What is missing is the
   genuinely-transformed subset. This is a real gap, not a complete one, and it
   has sat open across many ingestion cycles without blocking anything.
3. **The failure mode is already safe.** Such files produce
   `DicomError::CompressedPixelData` with an actionable message, not a wrong
   image. The cost of the gap is that a rare object needs decompressing at the
   archive boundary — a documented, supported workflow.

If Part 2 support is later judged necessary, the decision order should be:

1. **Obtain ISO/IEC 15444-2.** Everything else is gated on this. If the
   organisation already licenses it, that removes the largest obstacle
   immediately.
2. **Scope it as a fork or an upstream contribution to `pdfluent-jpeg2000`,
   not a second codec.** The existing crate already has the MQ coder, tier-2
   and Part 1 MCT that Part 2 builds on; the delta is marker parsing plus the
   transform. A parallel decoder here would mean maintaining a second JPEG 2000
   implementation.
3. **Only then consider OpenJPEG**, and treat it as a workspace policy change
   (C dependency, WASM port) rather than a crate detail.

## What this spike did not settle

- Whether Part 2 multi-component is common enough in the image population this
  project actually serves to justify the work. That is a question about real
  archives, not about the standard, and it should be answered with data before
  anyone funds the standard's purchase.
- Whether `pdfluent-jpeg2000`'s maintainers would accept a Part 2 contribution.
  Cheap to ask, and worth asking early — it is the cheapest of the three routes
  to rule in or out.
