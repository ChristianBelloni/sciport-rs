pub fn main() {
    #[cfg(all(target_os = "macos", feature = "blas"))]
    println!("cargo:rustc-link-lib=framework=Accelerate");
}
