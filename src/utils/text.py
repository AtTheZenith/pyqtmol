import re


def subscript_formula(formula: str) -> str:
    """
    Converts a chemical formula string to use subscript numbers.
    Example: H2O -> H₂O
    """
    subscript_map = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    match = re.match(r"^(\d+)?([A-Za-z].*)$", formula)
    if not match:
        return formula

    coeff, compound = match.groups()

    parts = []
    i = 0
    while i < len(compound):
        c = compound[i]
        if c.isalpha() or c in "()":
            parts.append(c)
            i += 1
        elif c.isdigit():
            j = i
            while j < len(compound) and compound[j].isdigit():
                j += 1
            subscript = compound[i:j].translate(subscript_map)
            parts.append(subscript)
            i = j
        else:
            parts.append(c)
            i += 1

    transformed = "".join(parts)
    return (coeff or "") + transformed
