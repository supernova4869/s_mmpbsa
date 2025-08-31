use xdrfile::*;
use indicatif::{ProgressBar, ProgressStyle};

pub fn read_xtc(xtc: &str) -> Vec<(f32, Vec<[f32; 3]>)> {
    // 读取所有帧，但只提取需要的数据
    let trj = XTCTrajectory::open_read(xtc).unwrap();
    
    // 先收集所有坐标数据的引用，避免持有整个 Frame
    let read_pb = ProgressBar::new_spinner();
    read_pb.set_style(
        ProgressStyle::default_spinner()
            .tick_chars("⠁⠂⠄⡀⢀⠠⠐⠈ ")
            .template("{spinner} Reading trajectory: {msg} {elapsed_precise}").unwrap()
    );
    let mut frame_count = 0;

    let frame_data: Vec<_> = trj.into_iter()
        .map(|result| {
            frame_count += 1;
            read_pb.inc(1);
            read_pb.set_message(format!("Read {} frames", frame_count));
            result.map(|frame| (frame.time, frame.coords.to_vec()))
        })
        .collect::<Result<Vec<(f32, Vec<[f32; 3]>)>, _>>().unwrap();
    
    println!("Finished reading trajectory with {} frames.", frame_data.len());

    frame_data
}