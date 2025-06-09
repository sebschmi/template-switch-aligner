use resvg::usvg::Options;

pub mod error;
pub mod plain_text;
pub mod svg;
pub mod ts_arrangement;

pub fn svg_to_png(svg_in: &[u8], zoom: f32) -> Vec<u8> {
    log::info!("Converting SVG to PNG");

    let mut options = Options::default();
    options.fontdb_mut().load_system_fonts();
    let svg = resvg::usvg::Tree::from_data(svg_in, &options).unwrap();

    let raster_image_size = svg.size();
    let raster_image_width = (raster_image_size.width().ceil() * zoom) as u32;
    let raster_image_height = (raster_image_size.height().ceil() * zoom) as u32;
    log::info!("PNG size: {raster_image_width}x{raster_image_height}",);

    let mut raster_image =
        resvg::tiny_skia::Pixmap::new(raster_image_width, raster_image_height).unwrap();
    resvg::render(
        &svg,
        resvg::usvg::Transform::from_scale(zoom, zoom),
        &mut raster_image.as_mut(),
    );
    raster_image.encode_png().unwrap()
}
