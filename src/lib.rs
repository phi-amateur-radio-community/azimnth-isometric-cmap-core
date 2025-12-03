use shapefile::{Shape, ShapeReader};
use std::io::Cursor;
use std::mem::ManuallyDrop;
use wide::{f64x4, f64x8};
use rayon::prelude::*;
use std::iter::Zip;
use wide::CmpLt;

const DEGREE_TO_RADIAN_CONSTANT: f64 = std::f64::consts::PI / 180_f64;
const RADIAN_TO_DEGREE_CONSTANT: f64 = 180_f64 / std::f64::consts::PI;
const FRAC_PI_DEGREE: f64x8 = f64x8::new([180f64; 8]); 

// TODO C support of ltargcm
//pub fn latlon_to_azimnth_isometric_csupport() {}
pub fn latlon_to_azimnth_isometric(
    latitude_delta: f64,
    longitude_delta: f64,
) -> (f64, f64) {
    let square = |x: f64| x * x;
    let degree_to_radian = |degree: f64| degree * DEGREE_TO_RADIAN_CONSTANT;
    let radian_to_degree = |radian: f64| radian * RADIAN_TO_DEGREE_CONSTANT;
    let latitude_delta_radian: f64 = degree_to_radian(latitude_delta);
    let longitude_delta_radian: f64 = degree_to_radian(longitude_delta);
    let hemispheres_anterior_or_posterior: bool = longitude_delta.abs() < 90_f64;

    #[cfg(test)]
    dbg!(latitude_delta_radian, longitude_delta_radian);

    let latitude_delta_radian_sine: f64 = latitude_delta_radian.sin();
    let longitude_delta_radian_sine: f64 = longitude_delta_radian.sin();

    #[cfg(test)]
    dbg!(latitude_delta_radian_sine, longitude_delta_radian_sine);

    let distance: f64 =
        (square(latitude_delta_radian_sine) + square(longitude_delta_radian_sine)).sqrt();
    let k: f64 = match distance {
        0_f64 => 0_f64,
        _ => {
            (if hemispheres_anterior_or_posterior {
                radian_to_degree(distance.asin())
            } else {
                180_f64 - radian_to_degree(distance.asin())
            }) / distance
        }
    };

    #[cfg(test)]
    dbg!(distance, k);

    (
        longitude_delta_radian_sine * k,
        latitude_delta_radian_sine * k,
    )
}

/*
struct LatlonToAzimnthIsometricSimdIterator<'a> {
    latitude_delta_vec: Vec<f64>,
    longitude_delta_vec: Vec<f64>,
    chunks_bodys: Zip<ChunksExact<'a, f64>, ChunksExact<'a, f64>>,
    chunks_remainder: (&'a [f64], &'a [f64]),
}

impl LatlonToAzimnthIsometricSimdIterator<'_> {
    fn new(latitude_delta_vec: Vec<f64>, longitude_delta_vec: Vec<f64>) -> Self {
        let len = std::cmp::min(latitude_delta_vec.len(), longitude_delta_vec.len());
        let latitude_delta_vec_chunks = latitude_delta_vec[..len].chunks_exact(8);
        let longitude_delta_vec_chunks = longitude_delta_vec[..len].chunks_exact(8);
        LatlonToAzimnthIsometricSimdIterator {
            latitude_delta_vec,
            longitude_delta_vec,
            chunks_bodys: latitude_delta_vec_chunks.zip(longitude_delta_vec_chunks),
            chunks_remainder: (
                latitude_delta_vec_chunks.remainder(),
                longitude_delta_vec_chunks.remainder()
            ),
        }
    }

    fn remainder<'a>(&self) -> (&'a [f64], &'a [f64]) { self.chunks_remainder }
}

impl Iterator for LatlonToAzimnthIsometricSimdIterator<'_> {
    type Item = (Vec<f64>, Vec<f64>);

    fn next(&mut self) -> Option<Self::Item> {
        let square = |x: f64x8| x * x;
        let chunks_body = match self.chunks_bodys.next() {
            Some(s) => s,
            None => return None,
        };
        let chunks_latitude = f64x8::new(chunks_body.0.try_into().expect("Why the len of the chunks(8)'s item not 8?"));
        let chunks_longitude = f64x8::new(chunks_body.1.try_into().expect("Why the len of the chunks(8)'s item not 8?"));
        
        let chunks_posistion = chunks_longitude.abs() < f64x8::FRAC_PI_2;
        
        let chunks_latitude_radian = chunks_latitude.to_radians();
        let chunks_longitude_radian = chunks_longitude.to_radians();
        
        let chunks_latitude_radian_sine = chunks_latitude_radian.sin();
        let chunks_longitude_radian_sine = chunks_longitude_radian.sin();
        
        let distance = (square(chunks_latitude_radian_sine) + square(chunks_longitude_radian_sine)).sqrt();
        
        let distance_arcsine = distance.asin().to_degrees();
        let distance_arcsine_back = FRAC_PI_DEGREE - distance_arcsine;
        
        let k_uncheck = chunks_posistion.blend(distance_arcsine, distance_arcsine_back) / distance;

        let k_infinite_mask = k_uncheck.is_infinite();
        let k_nan_mask = k_uncheck.is_non();
        
        let k_bad = k_infinite_mask | k_nan_mask;
        
        let k = k_bad.blend(f64x8::ZERO, k_uncheck);

        Some((
//            (chunks_latitude * k).to_vec(),
//            (chunks_longitude * k).to_vec(),
        ))
    }
}
*/

// TODO used wide to accelerate multipoint,line,polygon
//pub fn latlon_to_azimnth_isometric_simd_csupport() {}

pub fn latlon_to_azimnth_isometric_simd(latitude_delta_vec: Vec<f64>, longitude_delta_vec: Vec<f64>) -> (Vec<f64>, Vec<f64>) {
    let len = std::cmp::min(latitude_delta_vec.len(), longitude_delta_vec.len());
    let latitude_delta_vec_chunks = latitude_delta_vec[..len].par_chunks_exact(8);
    let longitude_delta_vec_chunks = longitude_delta_vec[..len].par_chunks_exact(8);
    let chunks_remainder = latitude_delta_vec_chunks.remainder().iter().copied().zip(longitude_delta_vec_chunks.remainder().iter().copied());
    let chunks_bodys: Vec<([f64; 8], [f64; 8])> = latitude_delta_vec_chunks
        .zip(longitude_delta_vec_chunks)
        .map(| ( latitude_slice, longitude_slice ): (&[f64], &[f64]) | -> ([f64; 8], [f64; 8]) {
            let latitude_array = latitude_slice;
            let longitude_array = longitude_slice;
            
            let square = |x: f64x8| x * x;

            let latitude = f64x8::from(latitude_array);
            let longitude = f64x8::from(longitude_array);

            let posistion_array: [f64; 8] = longitude.abs().to_array();

//            for i in posistion_array {
//                posision_array

            let posistion = longitude.abs().simd_lt(f64x8::FRAC_PI_2);

            let latitude_radian = latitude.to_radians();
            let longitude_radian = longitude.to_radians();

            let latitude_radian_sine = latitude_radian.sin();
            let longitude_radian_sine = longitude_radian.sin();

            let distance = (square(latitude_radian_sine) + square(longitude_radian_sine)).sqrt();

            let distance_arcsine = distance.asin().to_degrees();
            let distance_arcsine_back = FRAC_PI_DEGREE - distance_arcsine;

            let k_uncheck = posistion.blend(distance_arcsine, distance_arcsine_back) / distance;

            let k_infinite_mask = k_uncheck.is_inf();
            let k_nan_mask = k_uncheck.is_nan();

            let k_bad = k_infinite_mask | k_nan_mask;

            let k = k_bad.blend(f64x8::ZERO, k_uncheck);

        (
            (latitude * k).to_array(),
            (longitude * k).to_array()
        )
    } ).collect();

    let mut xs: Vec<f64> = Vec::with_capacity(len);
    let mut ys: Vec<f64> = Vec::with_capacity(len);

    for (x, y) in chunks_bodys {
        xs.extend(x);
        ys.extend(y);
    }


    for (latitude, longitude) in chunks_remainder {
        let result = latlon_to_azimnth_isometric(latitude, longitude);
        xs.push(result.0);
        ys.push(result.1);
    }

    (xs, ys)
}

pub struct ReturnContent {
    pub status: bool,
    pub ptr: *const u8,
    pub len: usize,
}

impl ReturnContent {
    fn new(data: Vec<u8>, status: bool) -> Self {
        let data = ManuallyDrop::new(data);
        ReturnContent {
            status: status,
            ptr: data.as_ptr(),
            len: data.len(),
        }
    }
}

pub struct ColorData {
    pub color_point: u8,
    pub color_multipoint: u8,
    pub color_line: u8,
    pub color_polygon_line: u8,
    pub color_polygon_fill: u8,
}

pub fn shapefile_generate(
    buffer_ptr: *const u8,
    buffer_len: usize,
    color: ColorData,
) -> ReturnContent {
    let buffer = unsafe { std::slice::from_raw_parts(buffer_ptr, buffer_len) };
    let cursor = Cursor::new(buffer);
    let mut reader = match ShapeReader::new(cursor) {
        Ok(data) => data,
        Err(e) => return ReturnContent::new(e.to_string().into_bytes(), false),
    };
    // TODO Draw the picture used par_iter
    //    reader.for_each(|shape| )
    ReturnContent::new(Vec::new(), true)
}

// TODO Draw fuction
//fn shapefile_draw() -> {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn transprojection() {
        let result = latlon_to_azimnth_isometric(0_f64, 0_f64);
        assert_eq!(result, (0_f64, 0_f64));

        let result = latlon_to_azimnth_isometric(-90_f64, 0_f64);
        assert_eq!(result, (0_f64, -90_f64));

        let result = latlon_to_azimnth_isometric(0_f64, 180_f64);
        assert_eq!(result, (180_f64, 0_f64));
    }
}
