use std::collections::HashSet;

const TOPIC_LEN: usize = 20;
// don't modify across threads
static mut DEBUG_TOPICS: Option<HashSet<&str>> = None;

#[macro_export]
#[cfg(debug_assertions)]
macro_rules! debugp {
    ($topic:expr, $($arg:tt)*) => (
        if $crate::debug_print::is_topic_enabled($topic) {
            let topic_padded = $crate::debug_print::pad_topic($topic);
            println!("{}{}", topic_padded, format!($($arg)*))
        }
    );
}

#[macro_export]
#[cfg(not(debug_assertions))]
macro_rules! debugp {
    ($topic:expr, $($arg:tt)*) => {};
}

pub fn is_topic_enabled(topic: &str) -> bool {
    unsafe {
        match DEBUG_TOPICS {
            Some(ref debug_topics) => debug_topics.contains(topic),
            // print all debug by default to have debug enabled in tests
            None => true,
        }
    }
}

pub fn pad_topic(topic: &str) -> String {
    let mut topic_str = topic.to_owned();
    topic_str.push(':');
    topic_str.push(' ');

    while topic_str.len() < TOPIC_LEN - 1 {
        topic_str.push(' ');
    }

    return topic_str;
}

pub fn enable_topic(topic: &'static str) {
    unsafe {
        if DEBUG_TOPICS.is_none() {
            DEBUG_TOPICS = Some(HashSet::new());
        }

        DEBUG_TOPICS.as_mut().unwrap().insert(topic);
    }
}

pub fn disable_all_topics() {
    unsafe {
        DEBUG_TOPICS = Some(HashSet::new());
    }
}
