//! Wrap an `impl Extend` to map elements before extending.

use std::marker::PhantomData;

/// Wraps a type `Extender: Extend<ExtendElement>` and applies a mapping function before passing elements to the `Extender`.
pub struct ExtendMap<
    'extender,
    InputElement,
    ExtendElement,
    Mapper: FnMut(InputElement) -> ExtendElement,
    Extender: Extend<ExtendElement>,
> {
    extender: &'extender mut Extender,
    mapper: Mapper,
    phantom_data: PhantomData<(InputElement, ExtendElement)>,
}

impl<
    'extender,
    InputElement,
    ExtendElement,
    Mapper: FnMut(InputElement) -> ExtendElement,
    Extender: Extend<ExtendElement>,
> ExtendMap<'extender, InputElement, ExtendElement, Mapper, Extender>
{
    pub fn new(extender: &'extender mut Extender, mapper: Mapper) -> Self {
        Self {
            extender,
            mapper,
            phantom_data: PhantomData,
        }
    }

    pub fn into_inner(self) -> &'extender mut Extender {
        self.extender
    }
}

impl<
    InputElement,
    ExtendElement,
    Mapper: FnMut(InputElement) -> ExtendElement,
    Extender: Extend<ExtendElement>,
> Extend<InputElement> for ExtendMap<'_, InputElement, ExtendElement, Mapper, Extender>
{
    fn extend<T: IntoIterator<Item = InputElement>>(&mut self, iter: T) {
        self.extender.extend(iter.into_iter().map(&mut self.mapper));
    }
}
