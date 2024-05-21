pub fn main() {
    #[cfg(all(target_os = "macos", feature = "blas"))]
    println!("cargo:rustc-link-lib=framework=Accelerate");
    #[cfg(target_os = "macos")]
    println!(
        "cargo:rustc-link-arg=-Wl,-rpath,/Library/Developer/CommandLineTools/Library/Frameworks"
    );
}
