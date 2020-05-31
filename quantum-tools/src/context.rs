use crate::system::System;

#[derive(Clone, Debug)]
pub struct Context<'a, T> {
    pub value: T,
    pub system: Option<&'a System>,
}
