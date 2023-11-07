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
    # print(a, b, c)

    # if choice == 1:
    #     (a, b, c), _ = op.curve_fit(
    #         func1, x, y
    #     )
    # elif choice == 2:
    #     (a, b, c), _ = op.curve_fit(
    #         func2, x, y
    #     )
    # elif choice == 3:
    #     (a, b, c), _ = op._curve_fit(
    #         func3, x, y
    #     )

    def glyc_update(i):
        new_glyc_value = functions[choice-1](i, a, b, c)
        return new_glyc_value * 100

    return glyc_update


if __name__ == "__main__":
    from import_elife import import_data

    light = import_data("light")
    light_calc = import_data("light calc")
    glycogen = (import_data("glycogen")/import_data("DW"))
    print(glycogen)
    fit(light, glycogen)
    fit(light_calc, glycogen)
