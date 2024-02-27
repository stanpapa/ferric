use crate::linear_algebra::{matrix::FMatrix, traits::Mat};

use std::collections::BTreeMap;
use std::ops::{Deref, DerefMut};
use std::{
    collections::btree_map::Entry,
    fmt::{Display, Formatter},
    fs::File,
    io::prelude::*,
    ops::{Index, IndexMut},
};

use serde::{de::DeserializeOwned, Deserialize, Deserializer, Serialize};

pub type FMatrixContainer = MatrixContainer<FMatrix>;

// #[derive(Clone, Deserialize, Serialize)]
#[derive(Clone)]
pub struct MatrixContainer<T: Mat> {
    data: BTreeMap<(usize, usize), T>,
}

impl<T: Mat> Default for MatrixContainer<T> {
    fn default() -> Self {
        Self {
            data: BTreeMap::default(),
        }
    }
}

impl<T: Mat> MatrixContainer<T> {
    pub fn new() -> MatrixContainer<T> {
        MatrixContainer::<T>::default()
    }

    pub fn insert(&mut self, index: (usize, usize), mat: &T) {
        self.data.insert(index, mat.clone());
    }

    // pub fn keys(&self) -> Keys<'_, (usize, usize), T> {
    //     self.data.keys()
    // }
    //
    // pub fn get(&self, key: &(usize, usize)) -> Option<&T> {
    //     self.data.get(key)
    // }

    pub fn get_mut(&mut self, key: &(usize, usize)) -> Option<&mut T> {
        self.data.get_mut(key)
    }

    pub fn entry(&mut self, key: &(usize, usize)) -> Entry<'_, (usize, usize), T> {
        self.data.entry(*key)
    }
}

impl<T: Mat> Deref for MatrixContainer<T> {
    type Target = BTreeMap<(usize, usize), T>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl<T: Mat> DerefMut for MatrixContainer<T> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

/// Functions for Serializing and Deserializing a MatrixContainer
#[derive(Serialize)]
struct MatrixContainerEntry<'entry, T: Mat + Serialize> {
    key: &'entry (usize, usize),
    val: &'entry T,
}
#[derive(Deserialize)]
struct MatrixContainerEntryOwned<T: Mat> {
    key: (usize, usize),
    val: T,
}

// impl<T: Mat + Serialize + DeserializeOwned> Into<Vec<MatrixContainerEntry<T>>>
//     for MatrixContainer<T>
// {
//     fn into(self) -> Vec<MatrixContainerEntry<T>> {
//         self.data
//             .iter()
//             .map(|(key, val)| MatrixContainerEntry {
//                 key: *key,
//                 val: val.clone(),
//             })
//             .collect()
//     }
// }

impl<T: Mat + Serialize> Serialize for MatrixContainer<T> {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        serializer.collect_seq(
            self.data
                .iter()
                .map(|(key, val)| MatrixContainerEntry { key, val }),
        )
    }
}

impl<'de, T: Mat + Deserialize<'de>> Deserialize<'de> for MatrixContainer<T> {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        Vec::<MatrixContainerEntryOwned<T>>::deserialize(deserializer).map(|mut v| Self {
            data: v.drain(..).map(|kv| (kv.key, kv.val)).collect(),
        })
    }
}

impl<T: Mat + Serialize + DeserializeOwned> MatrixContainer<T> {
    pub fn store(&self, name: &str) {
        let mut buffer = File::create(name).expect("Unable to create file");
        write!(
            buffer,
            "{}",
            serde_json::to_string(self).expect("Unable to serialize MatrixContainer")
        )
        .expect("Unable to write to file");
    }

    pub fn retrieve(name: &str) -> Self {
        let mut file = File::open(name).expect("Unable to open file for reading");
        let mut buffer = String::new();
        file.read_to_string(&mut buffer)
            .expect("Unable to read file");
        serde_json::from_str(&buffer).expect("Unable to deserialize a MatrixContainer")
    }
}

impl<T: Mat + Display> Display for MatrixContainer<T> {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        for (index, mat) in &self.data {
            writeln!(f, "{:?}\n{}", index, mat)?;
        }

        Ok(())
    }
}

impl<T: Mat> Index<(usize, usize)> for MatrixContainer<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        match self.data.get(&index) {
            Some(mat) => mat,
            None => panic!("No matrix found at {:?}", index),
        }
    }
}

impl<T: Mat> IndexMut<(usize, usize)> for MatrixContainer<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        match self.data.get_mut(&index) {
            Some(mat) => mat,
            None => panic!("No matrix found at {:?}", index),
        }
    }
}
