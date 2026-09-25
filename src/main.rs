// Hide the console window in release builds on Windows.
#![cfg_attr(all(not(debug_assertions), target_os = "windows"), windows_subsystem = "windows")]

mod db;
mod demo;
mod display;
mod fasta;
mod mapper;
mod model;
mod primer;
mod search;
mod seq;
mod thermo;
mod ui;

fn main() -> iced::Result {
    ui::run()
}
