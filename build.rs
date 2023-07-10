fn main() {
    println!(
        "cargo:rustc-link-search={}",
        "/Volumes/Babylon/SideProjects/sci-rs/libs/lib/"
    );
    println!("cargo:rustc-link-lib=dylib=complex_bessel");
}
