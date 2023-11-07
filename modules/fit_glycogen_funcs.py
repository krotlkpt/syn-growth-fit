from scipy import optimize as op


def fit(x, y, choice=2):
    def func1(input, a, b, c):
        result = a + input * b
        return result

    def func2(input, a, b, c):
        result = (a * input / (b + input) + c)
        return result

    def func3(input, a, b, c):
        result = b * -input + c
        return result

    functions = [func1, func2, func3]

    (a, b, c), _ = op.curve_fit(
        functions[choice-1],
        x, y,
        p0=(1, 170, 0.05)
    )

    def glyc_update(i):
        new_glyc_value = functions[choice-1](i, a, b, c)
        return new_glyc_value * 100

    return glyc_update
