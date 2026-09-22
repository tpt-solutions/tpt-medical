#[test]
fn debug_walk() {
    let series = crate::synthetic::femur_phantom(8, 8, 4);
    let bytes = &series.slices[0].bytes;
    let mut pos = 132;
    while pos + 8 <= bytes.len() {
        let g = u16::from_le_bytes([bytes[pos], bytes[pos+1]]);
        let e = u16::from_le_bytes([bytes[pos+2], bytes[pos+3]]);
        let vr = core::str::from_utf8(&bytes[pos+4..pos+6]).unwrap_or("??");
        let (len, hdr) = if matches!(vr, "OB"|"OW"|"SQ"|"UN"|"UT"|"OL"|"OV"|"UC"|"SV"|"UV") {
            (u32::from_le_bytes([bytes[pos+8],bytes[pos+9],bytes[pos+10],bytes[pos+11]]) as usize, 12)
        } else {
            (u16::from_le_bytes([bytes[pos+6],bytes[pos+7]]) as usize, 8)
        };
        println!("pos={pos} tag=({g:04x},{e:04x}) vr={vr} len={len}");
        pos += hdr + len;
        if pos > bytes.len() { println!("OVERRUN to {pos} (file {})", bytes.len()); break; }
    }
    println!("end pos={pos} file_len={}", bytes.len());
    panic!("debug done");
}
