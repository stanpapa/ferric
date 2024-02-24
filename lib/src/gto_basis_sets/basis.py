import json
import requests
import argparse

element_names = {
    "1": "Hydrogen",
    "2": "Helium",
    "3": "Lithium",
    "4": "Beryllium",
    "5": "Boron",
    "6": "Carbon",
    "7": "Nitrogen",
    "8": "Oxygen",
    "9": "Fluorine",
    "10": "Neon",
    "11": "Sodium",
    "12": "Magnesium",
    "13": "Aluminum",
    "14": "Silicon",
    "15": "Phosphorus",
    "16": "Sulfur",
    "17": "Chlorine",
    "18": "Argon",
    "19": "Potassium",
    "20": "Calcium",
    "21": "Scandium",
    "22": "Titanium",
    "23": "Vanadium",
    "24": "Chromium",
    "25": "Manganese",
    "26": "Iron",
    "27": "Cobalt",
    "28": "Nickel",
    "29": "Copper",
    "30": "Zinc",
    "31": "Gallium",
    "32": "Germanium",
    "33": "Arsenic",
    "34": "Selenium",
    "35": "Bromine",
    "36": "Krypton",
    "37": "Rubidium",
    "38": "Strontium",
    "39": "Yttrium",
    "40": "Zirconium",
    "41": "Niobium",
    "42": "Molybdenum",
    "43": "Technetium",
    "44": "Ruthenium",
    "45": "Rhodium",
    "46": "Palladium",
    "47": "Silver",
    "48": "Cadmium",
    "49": "Indium",
    "50": "Tin",
    "51": "Antimony",
    "52": "Tellurium",
    "53": "Iodine",
    "54": "Xenon",
    "55": "Cesium",
    "56": "Barium",
    "57": "Lanthanum",
    "58": "Cerium",
    "59": "Praseodymium",
    "60": "Neodymium",
    "61": "Promethium",
    "62": "Samarium",
    "63": "Europium",
    "64": "Gadolinium",
    "65": "Terbium",
    "66": "Dysprosium",
    "67": "Holmium",
    "68": "Erbium",
    "69": "Thulium",
    "70": "Ytterbium",
    "71": "Lutetium",
    "72": "Hafnium",
    "73": "Tantalum",
    "74": "Tungsten",
    "75": "Rhenium",
    "76": "Osmium",
    "77": "Iridium",
    "78": "Platinum",
    "79": "Gold",
    "80": "Mercury",
    "81": "Thallium",
    "82": "Lead",
    "83": "Bismuth",
    "84": "Polonium",
    "85": "Astatine",
    "86": "Radon",
    "87": "Francium",
    "88": "Radium",
    "89": "Actinium",
    "90": "Thorium",
    "91": "Protactinium",
    "92": "Uranium",
    "93": "Neptunium",
    "94": "Plutonium",
    "95": "Americium",
    "96": "Curium",
    "97": "Berkelium",
    "98": "Californium",
    "99": "Einsteinium",
    "100": "Fermium",
    "101": "Mendelevium",
    "102": "Nobelium",
    "103": "Lawrencium",
    "104": "Rutherfordium",
    "105": "Dubnium",
    "106": "Seaborgium",
    "107": "Bohrium",
    "108": "Hassium",
    "109": "Meitnerium",
    "110": "Darmstadtium",
    "111": "Roentgenium",
    "112": "Copernicium",
    "113": "Nihonium",
    "114": "Flerovium",
    "115": "Moscovium",
    "116": "Livermorium",
    "117": "Tennessine",
    "118": "Oganesson",
}

labels = {
    0: "s",
    1: "p",
    2: "d",
    3: "f",
}


# parse URL
def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert JSON from BasisSetExchange to Ferric source code"
    )
    parser.add_argument("url", type=str, help="URL to JSON")
    return parser.parse_args()


# retrieve JSON file from URL
def read_basis_set(url):
    try:
        response = requests.get(url)
        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            return response.text
        else:
            print(f"Error: {response.status_code}")
            return None
    except Exception as e:
        print(f"An error occurred: {e}")
        return None


# generate filtered list of exponents and coefficients
def exp_coef(shell, s):
    exps = [float(i) for i in shell["exponents"]]
    coeffs = [float(i) for i in shell["coefficients"][s]]
    return [e for e, c in zip(exps, coeffs) if c != 0.0], [
        c for c in coeffs if c != 0.0
    ]


def run():
    # read basis set from URL
    args = parse_args()
    json_raw = read_basis_set(args.url)

    # construct proper JSON object
    parsed_json = json.loads(json_raw)

    # basis set name
    name = parsed_json["name"]
    lowercase = name.lower().replace("-", "_")

    # sorted list of element numbers
    keys = [int(i) for i in parsed_json["elements"]]
    elements = [str(i) for i in sorted(keys)]

    # header
    print(
        f"""use crate::geometry::atom::Atom;
use crate::gto_basis_sets::basis::{{Basis, Shell}};

pub fn load_{lowercase}(atoms: &[Atom]) -> Basis {{
    println!("Loading {name} basis set");
    const MAX_ATOMIC_NUMBER: usize = {keys[-1]};

    let mut shells: Vec<Vec<Shell>> = vec![Default::default(); MAX_ATOMIC_NUMBER + 1];
"""
    )

    for el in elements:
        print(
            f"""
    // ----------------------------
    // Element #{el}, {element_names[el]}
    // ----------------------------
"""
        )
        shells = parsed_json["elements"][el]["electron_shells"]
        nshells = [len(shell["coefficients"]) for shell in shells]
        print(f"    shells[{el}] = Vec::with_capacity({sum(nshells)});\n")

        # initialise shells
        for shell in shells:
            nl = len(shell["angular_momentum"])
            nshell = len(shell["coefficients"])
            # 1 angular momentum, for every coefficient
            if nl == nshell:
                for s in range(0, nl):
                    exps, coeffs = exp_coef(shell, s)
                    print(
                        f"""    // {labels[shell["angular_momentum"][s]]}
    shells[{el}].push(Shell::new(
        {shell["angular_momentum"][s]},
        vec!{exps},
        vec!{coeffs},
    ));
"""
                    )
            # 1 angular momentum, multiple coefficients
            elif nl == 1 and nshell > 1:
                for s in range(0, nshell):
                    exps, coeffs = exp_coef(shell, s)
                    print(
                        f"""    // {labels[shell["angular_momentum"][0]]}
    shells[{el}].push(Shell::new(
        {shell["angular_momentum"][0]},
        vec!{exps},
        vec!{coeffs},
    ));
"""
                    )
            else:
                print("Unsupported shell layout")
                exit(1)

    print("    Basis::new(atoms, shells)\n}")


if __name__ == "__main__":
    run()
