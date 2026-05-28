from cdasws import CdasWs

dataset = "OMNI2_H0_MRG1HR"

cdas = CdasWs()

print("Inspecting OMNI dataset")
print(f"Dataset: {dataset}")
print()

datasets = cdas.get_datasets(id=dataset)

if not datasets:
    raise RuntimeError("Dataset was not found in CDAWeb.")

print("Dataset found.")
print()

variables = cdas.get_variables(dataset)

keywords = [
    "speed",
    "flow",
    "velocity",
    "density",
    "proton",
    "magnetic",
    "field",
    "imf",
    "bx",
    "by",
    "bz",
    "gsm",
    "gse",
    "magnitude",
]

print("Relevant variable candidates:")
print()

count = 0

for var in variables:
    text = " ".join(str(value) for value in var.values()).lower()

    if any(word in text for word in keywords):
        name = var.get("Name", "")
        short = var.get("ShortDescription", "")
        long = var.get("LongDescription", "")

        print(f"Name: {name}")
        if short:
            print(f"  Short: {short}")
        if long:
            print(f"  Long: {long[:250]}")
        print()

        count += 1

print(f"Printed {count} candidate variables.")
print()
print("Status: pass")
