from cdasws import CdasWs

cdas = CdasWs()

dataset = "OMNI2_H0_MRG1HR"

variables = [
    "V1800",
    "BX_GSE1800",
    "BY_GSE1800",
    "BZ_GSE1800",
    "ABS_B1800",
    "N1800",
]

start = "2024-07-25T00:00:00Z"
end = "2024-07-27T00:00:00Z"

status, data = cdas.get_data(dataset, variables, start, end)

print("status type:", type(status))
print("status:")
print(status)
print()

print("data type:", type(data))
print("data keys:")
print(list(data.keys()))
print()

for key, value in data.items():
    print("KEY:", key)
    print("  type:", type(value))

    if hasattr(value, "shape"):
        print("  shape:", value.shape)

    try:
        print("  len:", len(value))
    except Exception:
        pass

    text = repr(value)
    print("  repr first 500:")
    print(text[:500])
    print()
