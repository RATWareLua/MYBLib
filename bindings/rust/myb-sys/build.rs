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
    let directories: Vec<PathBuf> = env::var_os("MYB_LIB_DIR")
        .map(|value| env::split_paths(&value).collect())
        .unwrap_or_default();
    // dependents first: a static library is searched for what was asked for before it
    let mut libraries = Vec::new();
    for (feature, name) in [
        ("CARGO_FEATURE_SGP4", "myb_sgp4"),
        ("CARGO_FEATURE_NRLMSISE00", "myb_nrlmsise00"),
        ("CARGO_FEATURE_TROPOSPHERE", "myb_troposphere"),
        ("CARGO_FEATURE_IONOSPHERE", "myb_ionosphere"),
    ] {
        if env::var_os(feature).is_some() {
            libraries.push(name);
        }
    }
    libraries.push("myb");

    if directories.is_empty() {
        // an installed library: pkg-config knows where it is and what it needs
        let asked = Command::new("pkg-config").args(["--libs", "--static", "myb"]).output();
        if let Ok(output) = asked {
            if output.status.success() {
                for word in String::from_utf8_lossy(&output.stdout).split_whitespace() {
                    if let Some(directory) = word.strip_prefix("-L") {
                        println!("cargo:rustc-link-search=native={directory}");
                    }
                }
            }
        }
    }
    for directory in &directories {
        println!("cargo:rustc-link-search=native={}", directory.display());
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
    for name in &libraries {
        // the MSVC linker does not look for lib<name>.dll.a: name the import library itself
        let import = directories
            .iter()
            .map(|d| d.join(format!("lib{name}.dll.a")))
            .find(|path| path.exists());
        match import {
            Some(path) if kind == "dylib" && target_env == "msvc" => {
                println!("cargo:rustc-link-arg={}", path.display());
            }
            _ => println!("cargo:rustc-link-lib={kind}={name}"),
        }
    }
    if kind == "static" && target_os != "windows" {
        // what the static library itself needs: its own square roots, its threads
        println!("cargo:rustc-link-lib=m");
        println!("cargo:rustc-link-lib=pthread");
    }
}
