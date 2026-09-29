//! Ed25519 asymmetric signatures (RFC 8032), behind the `ed25519` cargo
//! feature — the **non-repudiation** counterpart to the symmetric HMAC
//! path: a signature made with a private key verifies against the public
//! key alone, so a signer cannot later disown a record by claiming a
//! shared key leaked.
//!
//! The primitive wraps `ed25519-dalek` (a reviewed implementation — this
//! crate deliberately does not hand-roll asymmetric crypto the way it
//! hand-rolls SHA-256/HMAC, which have compact, vector-verifiable
//! definitions). Keys are 32-byte seeds / compressed points; signatures
//! are 64 bytes. Verify against the public key only.

use ed25519_dalek::{Signature, Signer, Verifier};

/// An Ed25519 signing key (the 32-byte seed).
#[derive(Clone)]
pub struct SigningKey {
    inner: ed25519_dalek::SigningKey,
}

impl core::fmt::Debug for SigningKey {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        // The seed is the secret: never render it.
        f.debug_struct("SigningKey").finish_non_exhaustive()
    }
}

impl SigningKey {
    /// Derives the signing key from a 32-byte seed. Key generation (where
    /// the seed comes from) is the caller's concern — use a vetted CSPRNG;
    /// this crate does not generate keys.
    pub fn from_seed(seed: &[u8; 32]) -> Self {
        Self {
            inner: ed25519_dalek::SigningKey::from_bytes(seed),
        }
    }

    /// The public verification key to publish alongside records.
    pub fn verifying_key(&self) -> VerifyingKey {
        VerifyingKey {
            inner: self.inner.verifying_key(),
        }
    }

    /// Signs a message (e.g. a chain digest or a detached export tag),
    /// returning the 64-byte detached signature.
    pub fn sign(&self, message: &[u8]) -> [u8; 64] {
        self.inner.sign(message).to_bytes()
    }
}

/// An Ed25519 verification key (the 32-byte compressed point).
#[derive(Clone, Debug)]
pub struct VerifyingKey {
    inner: ed25519_dalek::VerifyingKey,
}

impl VerifyingKey {
    /// Parses a published public key.
    ///
    /// # Errors
    ///
    /// Returns the dalek error for a non-canonical/unreducible point
    /// encoding.
    pub fn from_bytes(bytes: &[u8; 32]) -> Result<Self, ed25519_dalek::SignatureError> {
        Ok(Self {
            inner: ed25519_dalek::VerifyingKey::from_bytes(bytes)?,
        })
    }

    /// The compressed public-key bytes.
    pub fn to_bytes(&self) -> [u8; 32] {
        self.inner.to_bytes()
    }

    /// Verifies a detached signature over a message. False for any
    /// malformed signature — never an error and never a panic.
    pub fn verify(&self, message: &[u8], signature: &[u8; 64]) -> bool {
        let Ok(sig) = Signature::from_slice(signature) else {
            return false;
        };
        self.inner.verify(message, &sig).is_ok()
    }
}

/// Convenience: verifies a chain digest under a published key — the shape
/// a forensic tool uses to check who signed an audit chain.
pub fn verify_chain_signature(
    key: &VerifyingKey,
    chain_digest_hex: &str,
    signature: &[u8; 64],
) -> bool {
    key.verify(chain_digest_hex.as_bytes(), signature)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{hash_chain, hex};

    /// RFC 8032 §7.1 TEST 2 (Ed25519): seed → public key derivation and
    /// the deterministic signature over the single-byte message `0x72` are
    /// checked against the standard's constants (cross-checked against an
    /// independent OpenSSL derivation when the constants were recorded).
    #[test]
    fn rfc8032_test_vector_2() {
        let seed = [
            0x4c, 0xcd, 0x08, 0x9b, 0x28, 0xff, 0x96, 0xda, 0x9d, 0xb6, 0xc3, 0x46, 0xec, 0x11,
            0x4e, 0x0f, 0x5b, 0x8a, 0x31, 0x9f, 0x35, 0xab, 0xa6, 0x24, 0xda, 0x8c, 0xf6, 0xed,
            0x4f, 0xb8, 0xa6, 0xfb,
        ];
        let key = SigningKey::from_seed(&seed);
        let vk_bytes = key.verifying_key().to_bytes();
        let expected_pk = [
            0x3d, 0x40, 0x17, 0xc3, 0xe8, 0x43, 0x89, 0x5a, 0x92, 0xb7, 0x0a, 0xa7, 0x4d, 0x1b,
            0x7e, 0xbc, 0x9c, 0x98, 0x2c, 0xcf, 0x2e, 0xc4, 0x96, 0x8c, 0xc0, 0xcd, 0x55, 0xf1,
            0x2a, 0xf4, 0x66, 0x0c,
        ];
        assert_eq!(vk_bytes, expected_pk, "public key derivation");
        let message = [0x72u8];
        let sig = key.sign(&message);
        let expected_sig = [
            0x92, 0xa0, 0x09, 0xa9, 0xf0, 0xd4, 0xca, 0xb8, 0x72, 0x0e, 0x82, 0x0b, 0x5f, 0x64,
            0x25, 0x40, 0xa2, 0xb2, 0x7b, 0x54, 0x16, 0x50, 0x3f, 0x8f, 0xb3, 0x76, 0x22, 0x23,
            0xeb, 0xdb, 0x69, 0xda, 0x08, 0x5a, 0xc1, 0xe4, 0x3e, 0x15, 0x99, 0x6e, 0x45, 0x8f,
            0x36, 0x13, 0xd0, 0xf1, 0x1d, 0x8c, 0x38, 0x7b, 0x2e, 0xae, 0xb4, 0x30, 0x2a, 0xee,
            0xb0, 0x0d, 0x29, 0x16, 0x12, 0xbb, 0x0c, 0x00,
        ];
        assert_eq!(sig, expected_sig, "deterministic signature");
        // Verification accepts the true signature and rejects tampering.
        let parsed = VerifyingKey::from_bytes(&expected_pk).expect("canonical point");
        assert!(parsed.verify(&message, &sig));
        assert!(!parsed.verify(b"not the message", &sig));
        let mut flipped = sig;
        flipped[0] ^= 0x01;
        assert!(!parsed.verify(&message, &flipped));
    }

    #[test]
    fn chain_digest_signing_round_trips() {
        // The non-repudiation workflow: sign the chain's head digest with
        // Ed25519, publish key + signature, verify without the secret.
        let entries = vec!["create:s1".to_string(), "export:s1".to_string()];
        let chain = hash_chain(b"run-ed", &entries);
        let head_hex = chain.last().expect("non-empty");
        let key = SigningKey::from_seed(&[42u8; 32]);
        let sig = key.sign(head_hex.as_bytes());
        let published = key.verifying_key();
        assert!(verify_chain_signature(&published, head_hex, &sig));
        // A different chain head fails.
        let other = hash_chain(b"run-other", &entries);
        assert!(!verify_chain_signature(
            &published,
            other.last().expect("non-empty"),
            &sig
        ));
        // Digests remain the plain SHA-256 chain (hex length sanity).
        assert_eq!(hex(&[0u8; 32]).len(), 64);
    }
}
