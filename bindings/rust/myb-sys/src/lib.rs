// SPDX-License-Identifier: BSD-3-Clause
//! The C interface of myb — every function, structure and constant of its public headers —
//! as Rust declarations, by the names of the headers. The declarations are generated from
//! the headers (`bindings/rust/tools/generate`, bindgen) and kept in the crate, so building
//! it needs neither libclang nor another crate.
//!
//! Enumerations are integers with constants (`myb_status`, `MYB_OK`, `MYB_E_DOMAIN`, ...).
//! Everything is `unsafe` to call, as C is; the crate `myb` wraps the path most programs
//! take in safe types.
//!
//! The optional libraries of `ext/` are modules behind features of the same names: `sgp4`,
//! `nrlmsise00`, `troposphere`, `ionosphere`, `jb2008`.
#![allow(non_camel_case_types, non_upper_case_globals, non_snake_case)]
#![allow(clippy::all)]

include!("bindings.rs");

#[cfg(feature = "sgp4")]
pub mod sgp4 {
    include!("ext_sgp4.rs");
}

#[cfg(feature = "nrlmsise00")]
pub mod nrlmsise00 {
    include!("ext_nrlmsise00.rs");
}

#[cfg(feature = "troposphere")]
pub mod troposphere {
    include!("ext_troposphere.rs");
}

#[cfg(feature = "ionosphere")]
pub mod ionosphere {
    include!("ext_ionosphere.rs");
}

#[cfg(feature = "jb2008")]
pub mod jb2008 {
    include!("ext_jb2008.rs");
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::ffi::CStr;

    #[test]
    fn the_library_is_the_one_declared() {
        // SAFETY: functions without preconditions; the name is a static string of the library
        unsafe {
            assert_eq!(myb_version_number() as u32, MYB_VERSION_NUMBER);
            let name = CStr::from_ptr(myb_status_name(MYB_E_FORMAT));
            assert_eq!(name.to_str().unwrap(), "MYB_E_FORMAT");
            let t = myb_time_from_seconds(12.5);
            assert_eq!((t.s, t.f), (12, 0.5));
        }
    }

    #[cfg(feature = "sgp4")]
    #[test]
    fn an_optional_library_is_linked() {
        let line1 = c"1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753";
        let line2 = c"2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";
        let mut tle = sgp4::myb_tle::default();
        let mut diag = myb_diag::default();
        // SAFETY: two NUL-terminated lines of 69 characters; valid places for the outputs
        let code = unsafe { sgp4::myb_tle_parse(line1.as_ptr(), line2.as_ptr(), 0, &mut tle, &mut diag) };
        assert_eq!(code, MYB_OK);
        assert_eq!(tle.catalog_number, 5);
        assert!((tle.inclination - 34.2682).abs() < 1e-12);
    }
}
