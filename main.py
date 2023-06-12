import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate


def ak(k, t):
    s = integrate.quad(lambda x: (3*x - 1) * np.cos(np.pi * k * x / 3), 0, 1)[0] * (1 / 3) +\
        integrate.quad(lambda x: 1 * np.cos(np.pi * k * x / 3), 1, 3)[0] * (1 / 3)
    if t == 1:
        s += integrate.quad(lambda x: f3(x) * np.cos(np.pi * k * x / 3), -3, -2)[0] * (1 / 3) +\
             integrate.quad(lambda x: f2(x) * np.cos(np.pi * k * x / 3), -2, 0)[0] * (1 / 3)
    else:
        s += integrate.quad(lambda x: (-f2(x) if t == 2 else f2(x)) * np.cos(np.pi * k * x / 3), -3, -1)[0] * (1 / 3) +\
             integrate.quad(lambda x: (f1(x) if t == 2 else f4(x)) * np.cos(np.pi * k * x / 3), -1, 0)[0] * (1 / 3)
    return s


def bk(k, t):
    s = integrate.quad(lambda x: (3 * x - 1) * np.sin(np.pi * k * x / 3), 0, 1)[0] * (1 / 3) + \
        integrate.quad(lambda x: 1 * np.sin(np.pi * k * x / 3), 1, 3)[0] * (1 / 3)
    if t == 1:
        s += integrate.quad(lambda x: f3(x) * np.sin(np.pi * k * x / 3), -3, -2)[0] * (1 / 3) + \
             integrate.quad(lambda x: f2(x) * np.sin(np.pi * k * x / 3), -2, 0)[0] * (1 / 3)
    else:
        s += integrate.quad(lambda x: (-f2(x) if t == 2 else f2(x)) * np.sin(np.pi * k * x / 3), -3, -1)[0] * (1 / 3) + \
             integrate.quad(lambda x: (f1(x) if t == 2 else f4(x)) * np.sin(np.pi * k * x / 3), -1, 0)[0] * (1 / 3)
    return s


def a0(t):
    s = integrate.quad(lambda x: 3 * x - 1, 0, 1)[0] * (1 / 3) + \
        integrate.quad(lambda x: 1, 1, 3)[0] * (1 / 3)
    if t == 1:
        s += integrate.quad(lambda x: f3(x), -3, -2)[0] * (1 / 3) + \
             integrate.quad(lambda x: f2(x), -2, 0)[0] * (1 / 3)
    else:
        s += integrate.quad(lambda x: -f2(x) if t == 2 else f2(x), -3, -1)[0] * (1 / 3) + \
             integrate.quad(lambda x: f1(x) if t == 2 else f4(x), -1, 0)[0] * (1 / 3)
    return s


def sigma(x, n, t):
    s = a0(t) / 2
    for k in range(1, n + 1):
        s += ak(k, t) * np.cos(np.pi * k * x / 3) + bk(k, t) * np.sin(np.pi * k * x / 3)
    return s


def f1(x):
    return 3 * x - 1


def f2(x):
    return 1 + 0 * x


def f3(x):
    return 3 * x + 8


def f4(x):
    return -3 * x - 1


def draw(axis, x1, x2, x3, x4, g1, g2, g3, g4):
    axis.plot(x1, g1)
    axis.plot(x2, g2)
    axis.plot(x3, g3)
    axis.plot(x4, g4)


def main():
    # Input data section
    t = int(input("Общий/по синусам/косинусам(1, 2, 3)?: "))
    n = int(input("Количество слагаемых (1-100): "))

    # plots
    fig, ax = plt.subplots()
    x1 = np.arange(0, 1, 0.01)
    x2 = np.arange(1, 3, 0.01)
    if t == 1:
        x3 = np.arange(-3, -2, 0.01)
        x4 = np.arange(-2, 0, 0.01)
        draw(ax, x1, x2, x3, x4, f1(x1), f2(x2), f3(x3), f2(x4))
    else:
        x3 = np.arange(-3, -1, 0.01)
        x4 = np.arange(-1, 0, 0.01)
        draw(ax, x1, x2, x3, x4, f1(x1), f2(x2), -f2(x3) if t == 2 else f2(x3), f1(x4) if t == 2 else f4(x4))
    x = np.arange(-3, 3, 0.01)
    ax.plot(x, sigma(x, n, t))
    # print(a0(1), a0(2), a0(3))
    ax.set(xlabel="x", ylabel="y", title=f"График функции и ряд Фурье (n = {n}, " +
                                         f"{'общий' if t == 1 else ('по синусам'if t == 2 else 'по косинусам')})")
    # ax.set(xlabel="x", ylabel="y", title=f'График функции, продолженной нечетным образом')
    ax.grid()
    plt.savefig(f"fig_{t}_{n}.png")
    # plt.savefig("fig.png")
    plt.show()


main()
