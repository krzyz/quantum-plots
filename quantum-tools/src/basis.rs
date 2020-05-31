use crate::context::Context;
use crate::system::System;

#[derive(Clone, Debug)]
pub enum Basis<'a, T> {
    Problem(Context<'a, T>),
    Solution(Context<'a, T>),
}

impl<T> Basis<'static, T> {
    pub fn with_system<'a>(self, system: &'a System) -> Basis<'a, T> {
        match self {
            Basis::Problem(context) => Basis::Problem(Context {
                value: context.value,
                system: Some(system),
            }),
            Basis::Solution(context) => Basis::Solution(Context {
                value: context.value,
                system: Some(system),
            }),
        }
    }
}

impl<'a, T> Basis<'a, T> {
    pub fn without_system(self) -> Basis<'static, T> {
        match self {
            Basis::Problem(context) => Basis::Problem(Context {
                value: context.value,
                system: None,
            }),
            Basis::Solution(context) => Basis::Solution(Context {
                value: context.value,
                system: None,
            }),
        }
    }

    pub fn get_value(&self) -> &T {
        match self {
            Basis::Problem(ref context) => &context.value,
            Basis::Solution(ref context) => &context.value,
        }
    }

    pub fn get_system(&self) -> Option<&'a System> {
        match self {
            Basis::Problem(ref context) => context.system,
            Basis::Solution(ref context) => context.system,
        }
    }
}

impl<T> PartialEq for Basis<'_, T>
where
    T: PartialEq,
{
    fn eq(&self, other: &Self) -> bool {
        match std::mem::discriminant(self) == std::mem::discriminant(other) {
            true => self.get_value() == other.get_value(),
            false => false,
        }
    }
}
