// SPDX-License-Identifier: BSD-3-Clause
//! Finds the library to link.
//!
//!   MYB_LIB_DIR  the directory of libmyb.a or of the shared library (several: separated as
//!                PATH is; the libraries of ext/ are looked for there too). Without it the
//!                answer of `pkg-config --libs myb` is used, and failing that the linker's
//!                own search.
//!   MYB_LINK     "static" or "dylib"; by default static if libmyb.a is in MYB_LIB_DIR.
//!
//! With the MSVC toolchain only the shared library can be linked (through its import
//! library): the static one is built for the GNU ABI.
use std::env;
use std::path::PathBuf;
use std::process::Command;

fn main() {
    println!("cargo:rerun-if-env-changed=MYB_LIB_DIR");
    println!("cargo:rerun-if-env-changed=MYB_LINK");
    let target_os = env::var("CARGO_CFG_TARGET_OS").unwrap_or_default();
    let target_env = env::var("CARGO_CFG_TARGET_ENV").unwrap_or_default();
    let target_vendor = env::var("CARGO_CFG_TARGET_VENDOR").unwrap_or_default();
    let directories: Vec<PathBuf> = env::var_os("MYB_LIB_DIR")
        .map(|value| env::split_paths(&value).collect())
        .unwrap_or_default();
    let out_dir = PathBuf::from(env::var_os("OUT_DIR").expect("cargo sets OUT_DIR"));
    // dependents first: a static library is searched for what was asked for before it
    let mut libraries = Vec::new();
    for (feature, name) in [
        ("CARGO_FEATURE_SGP4", "myb_sgp4"),
        ("CARGO_FEATURE_NRLMSISE00", "myb_nrlmsise00"),
        ("CARGO_FEATURE_TROPOSPHERE", "myb_troposphere"),
        ("CARGO_FEATURE_IONOSPHERE", "myb_ionosphere"),
        ("CARGO_FEATURE_JB2008", "myb_jb2008"),
    ] {
        if env::var_os(feature).is_some() {
            libraries.push(name);
        }
    }
    libraries.push("myb");

    let mut search = directories.clone();
    if directories.is_empty() {
        // an installed library: pkg-config knows where it is and what it needs
        let asked = Command::new("pkg-config").args(["--libs", "--static", "myb"]).output();
        if let Ok(output) = asked {
            if output.status.success() {
                for word in String::from_utf8_lossy(&output.stdout).split_whitespace() {
                    if let Some(directory) = word.strip_prefix("-L") {
                        search.push(PathBuf::from(directory));
                    }
                }
            }
        }
    }
    let has_static = directories.iter().any(|d| d.join("libmyb.a").exists());
    let kind = match env::var("MYB_LINK").as_deref() {
        Ok("static") => "static",
        Ok("dylib") => "dylib",
        Ok(other) => panic!("MYB_LINK must be \"static\" or \"dylib\", not {other:?}"),
        Err(_) if has_static && target_env != "msvc" => "static",
        Err(_) => "dylib",
    };
    if kind == "static" && target_env == "msvc" {
        panic!("libmyb.a is built for the GNU ABI: link the shared library (MYB_LINK=dylib)");
    }
    // The linker of Apple takes lib<name>.dylib before lib<name>.a in the same directory, and
    // rustc has no -Bstatic for it: where a build tree holds both, a binary linked with -l
    // (this crate's own tests; the crates that depend on it take the archives bundled into its
    // rlib) would get the shared library, with no path to load it from. The archives are copied
    // into OUT_DIR, where no dylib is, and that directory is searched first.
    if kind == "static" && target_vendor == "apple" {
        for name in &libraries {
            let archive = search
                .iter()
                .map(|d| d.join(format!("lib{name}.a")))
                .find(|path| path.exists());
            if let Some(path) = archive {
                std::fs::copy(&path, out_dir.join(format!("lib{name}.a")))
                    .expect("the archive could not be copied into OUT_DIR");
                println!("cargo:rerun-if-changed={}", path.display());
            }
        }
        println!("cargo:rustc-link-search=native={}", out_dir.display());
    }
    for directory in &search {
        println!("cargo:rustc-link-search=native={}", directory.display());
    }
    // The MSVC linker looks for <name>.lib, not for the lib<name>.dll.a of the GNU toolchain:
    // the import library is copied under that name into OUT_DIR. (Naming the file itself with
    // rustc-link-arg would reach only this crate's own binaries — such arguments do not pass to
    // the crates that depend on it, and their tests would miss every symbol.)
    let mut copied = false;
    for name in &libraries {
        let import = directories
            .iter()
            .map(|d| d.join(format!("lib{name}.dll.a")))
            .find(|path| path.exists());
        if let Some(path) = import {
            if kind == "dylib" && target_env == "msvc" {
                std::fs::copy(&path, out_dir.join(format!("{name}.lib")))
                    .expect("the import library could not be copied into OUT_DIR");
                println!("cargo:rerun-if-changed={}", path.display());
                copied = true;
            }
        }
        println!("cargo:rustc-link-lib={kind}={name}");
    }
    if copied {
        println!("cargo:rustc-link-search=native={}", out_dir.display());
    }
    if kind == "static" && target_os != "windows" {
        // what the static library itself needs: its own square roots, its threads
        println!("cargo:rustc-link-lib=m");
        println!("cargo:rustc-link-lib=pthread");
    }
}
