// Hide the console window in release builds on Windows.
#![cfg_attr(all(not(debug_assertions), target_os = "windows"), windows_subsystem = "windows")]

mod db;
mod demo;
mod display;
mod fasta;
mod genbank;
mod linked;
mod mapper;
mod model;
/// All network access (NCBI).
mod net;
mod primer;
mod search;
mod seq;
mod thermo;
mod tools;
mod ui;

fn main() -> iced::Result {
    ui::run()
}
