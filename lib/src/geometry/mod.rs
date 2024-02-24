pub mod atom;
pub mod molecule;

#[derive(Clone, Debug, PartialEq, Default)]
pub enum Unit {
    #[default]
    Ångström,
    AtomicUnits,
    Bohr,
}
