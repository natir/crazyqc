//! Crazyqc compute some simple quality control information on reads dataset faster as possible.
//!
//! Idea is parse file in main thread and use rayon to parallelize analyse of record

/* mod declaration block */
pub mod cli;
pub mod error;
pub mod input;
