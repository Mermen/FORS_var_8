import math
import matplotlib.pyplot as plt


def fun_for_T_wtf(T1, T2):
    T = T1 * T2 * math.log(T2 / T1) / (T2 - T1)
    return T


def fun_for_H_K_delta_B1(p, R, m1, m2, T1, T2, r0, W1, W2, kT, g):
    T_wtf = fun_for_T_wtf(T1, T2)
    M = (m1 + m2) / 2
    mediumT = (T1 + T2) / 2
    dT = T2 - T1
    ro = p * M / mediumT / R / 1000
    etta = 267 * math.pow(10, -8) * math.sqrt(M * T_wtf) / math.pow(r0, 2) / W2
    print("etta = ", etta)
    D = 0.01314 * math.pow(T_wtf, 2) / math.sqrt(M * T_wtf) / p / math.pow(r0, 2) / W1
    print("D = ", D)
    delta = (5.0 * math.factorial(9) * math.pow(etta, 2) * math.pow(D, 2) / math.pow(ro, 2) / math.pow(g, 2) / math.pow(
        dT / mediumT, 2)) ** (1.0 / 6.0)
    print("delta = ", delta)
    B1 = 20.0 * delta
    print("B1 = ", B1)
    alphaT = (m2 - m1) / (m2 + m1) * kT
    print("alphaT = ", alphaT)
    Kc = math.pow(ro, 3) * math.pow(g, 2) * math.pow(delta, 7) * B1 * math.pow(dT / mediumT, 2) / math.factorial(
        9) / math.pow(etta, 2) / D
    print("Kc = ", Kc)
    Kd = ro * D * delta * B1
    print("Kd = ", Kd)
    K = Kc + Kd
    print("K = ", K)
    H = alphaT * math.pow(ro, 2) * g * math.pow(delta, 3) * B1 * math.pow(dT / mediumT, 2) / math.factorial(6) / etta
    print("H = ", H)
    return [K, H, delta, B1]


def F1(c):
    f1 = (2 * c - 1) * math.log(c / (1 - c))
    return f1


def q(R, Rf):
    return R / Rf


def R(C):
    return C / (1 - C)


def concentrationpw(Rf, H1, K1, z):
    return (Rf * math.exp(H1 * z / (2 * K1))) / (1 + Rf * math.exp(H1 * z / (2 * K1)))


def omegapw(P, C, H1, Cz):
    return (2 * P * (C - Cz)) / (H1 * Cz * (1 - Cz))


def func_for_Ideal_Cascade_concentration_omega(Cf, Cp, Cw, H1, K1, P, W):
    Rp = R(Cp)
    Rw = R(Cw)
    Rf = R(Cf)
    qp = q(Rp, Rf)
    qw = q(Rw, Rf)
    Lp = 2 * K1 / H1 * math.log(qp)
    Lw = 2 * K1 / H1 * math.log(qw)
    N = 1000
    dz = (Lp - Lw) / (N - 1)
    z = [Lw + i * dz for i in range(N)]
    concentration = [concentrationpw(Rf, H1, K1, z[i]) for i in range(N)]
    plt.plot(z, concentration)
    plt.title("Распределение концентрации для идеального каскада C(z)")
    plt.xlabel("Z, m")
    plt.ylabel("C(z)")
    plt.grid()
    plt.savefig("Распределение концентрации для идеального каскада C(z).png", format='png', dpi=300)
    plt.show()
    omega = [omegapw(P if (z[i] > 0) else -W, Cp if (z[i] > 0) else Cw, H1, concentration[i]) for i in range(N)]
    omega1 = 0.75 * max(omega)
    Lpr = fun_for_Lw_rectangle(K1, H1, Cp, Cf, P, omega1)
    Lwr = fun_for_Lw_rectangle(K1, H1, Cw, Cf, W, omega1)
    dzr = (Lpr - Lwr) / (N - 1)
    zr = [Lwr + i * dzr for i in range(N)]
    omegar = [omega1 for i in range(N)]
    omegar[0] = 0
    omegar[999] = 0
    plt.plot(z, omega)
    plt.title("Распределение числа колон для идеального каскада")
    plt.xlabel("Z, m")
    plt.ylabel("w(z)")
    plt.grid()
    plt.savefig("Число идеальных колонн.png", format='png', dpi=300)
    plt.show()
    plt.plot(z, omega, "blue", label="ИК")
    plt.plot(zr, omegar, "red", label="ПК")
    plt.xlabel("Z, m")
    plt.ylabel("w(z)")
    plt.title("Число колонн")
    plt.legend()
    plt.grid()
    plt.savefig("Число колонн.png", format='png', dpi=300)
    plt.show()
    return max(omega)


def V(C):
    return (2 * C - 1) * math.log(C / (1 - C))


def fun_for_Lw_rectangle(K1, H1, Cw, Cf, W, omega1):
    psi = W / (omega1 * H1)
    alpha = 1 / 2 * (1 + psi) + math.sqrt(1 / 4 * (1 + psi) ** 2 - Cw * psi)
    betta = 1 / 2 * (1 + psi) - math.sqrt(1 / 4 * (1 + psi) ** 2 - Cw * psi)
    Lw = K1 / (H1 * (betta - alpha)) * math.log((Cw - alpha) / (Cf - alpha) * (Cf - betta) / (Cw - betta))
    return Lw


def fun_for_rectangle_con(K1, H1, Cw, Cf, W, omega1, z):
    psi = W / (omega1 * H1)
    alpha = 1 / 2 * (1 + psi) + math.sqrt(1 / 4 * (1 + psi) ** 2 - Cw * psi)
    betta = 1 / 2 * (1 + psi) - math.sqrt(1 / 4 * (1 + psi) ** 2 - Cw * psi)
    C = (alpha - betta * (Cf - alpha) / (Cf - betta) * math.exp(H1 / K1 * (betta - alpha) * z)) / (
            1 - (Cf - alpha) / (Cf - betta) * math.exp(H1 / K1 * (betta - alpha) * z))
    return C


def fun_for_plot_con_rectangle(K1, H1, Cw, Cp, Cf, P, W, omega1):
    N = 1000
    Lpr = fun_for_Lw_rectangle(K1, H1, Cp, Cf, P, omega1)
    Lwr = fun_for_Lw_rectangle(K1, H1, Cw, Cf, W, omega1)
    dzr = (Lpr - Lwr) / (N - 1)
    zr = [Lwr + i * dzr for i in range(N)]
    concentration = [
        fun_for_rectangle_con(K1, H1, Cp if (zr[i] > 0) else Cw, Cf, P if (zr[i] > 0) else -W, omega1, zr[i]) for i in
        range(N)]
    plt.plot(zr, concentration)
    plt.title("Распределение концентрации\n для оптимального прямоугольного каскада C(z)")
    plt.xlabel("Z, m")
    plt.ylabel("C(z)")
    plt.grid()
    plt.savefig("Распределение концентрации для оптимального прямоугольного каскада C(z).png", format='png', dpi=300)
    plt.show()
    return 0


def fun_for_B_opt_k_opt(Cf, Cp, Cw, F, H1, K1, P, W, delta, omega1):
    N = 100
    koefmin = 14
    koefmax = 17
    dkoef = (koefmax - koefmin) / (N - 1)
    koef = [koefmin + dkoef * i for i in range(N)]
    kpd = []
    for i in range(N):
        Lp = fun_for_Lw_rectangle(K1 * koef[i] / 20, H1 * koef[i] / 20, Cp, Cf, P, omega1)
        Lw = fun_for_Lw_rectangle(K1 * koef[i] / 20, H1 * koef[i] / 20, Cw, Cf, -W, omega1)
        Lambda_rectanangle = (Lp - Lw) * omega1
        Lambda_ideal_r = fun_for_Lambda_ideal(Cf, Cp, Cw, F, (H1 * koef[i]) / 20, (K1 * koef[i]) / 20, P, W)
        kpd.append(Lambda_ideal_r / Lambda_rectanangle)
    plt.plot(koef, kpd)
    plt.title("КПД формы")
    plt.xlabel("i")
    plt.ylabel("КПД(i)")
    plt.grid()
    plt.savefig("КПД формы.png", format='png', dpi=300)
    plt.show()
    print("kpd_max = ",max(kpd))
    k_opt = koef[kpd.index(max(kpd))]
    print("k_opt = ", k_opt)
    B_opt = k_opt * delta
    print("B_opt = ", B_opt)
    return 0


def fun_for_Lambda_ideal(Cf, Cp, Cw, F, H1, K1, P, W):
    Lambda_ideal = (P * V(Cp) + W * V(Cw) - F * V(Cf)) / (H1 ** 2 / (4 * K1))
    return Lambda_ideal


def main():
    p = 101325
    print("p = ", p)
    R = 8.314
    print("R = ", R)
    g = 9.81
    print("g = ", g)
    Cp = 0.9
    print("Cp = ", Cp)
    Cf = 0.0037
    print("Cf = ", Cf)
    Cw = 0.5 * Cf
    print("Cw = ", Cw)
    E = 3.8 * 100
    print("E = ", E)
    m1 = 28  # N14-N14
    print("m1 = ", m1)
    m2 = 29  # N14-N15
    print("m2 = ", m2)
    M = (m1 + m2) / 2
    print("M = ", M)
    T1 = 300
    print("T1 = ", T1)
    T2 = 600
    print("T2 = ", T2)
    mediumT = (T1 + T2) / 2
    print("mediumT = ", mediumT)
    ro = p * M / mediumT / R / 1000
    print("ro = ", ro)
    P = 1 * ro / 86400 / 1000
    print("P = ", P)
    Tx = 91.5
    print("Tx = ", Tx)
    r0 = 3.68
    print("r0 = ", r0)
    T = fun_for_T_wtf(T1, T2)
    print("T = ", T)
    print("T/Tx = ", T / Tx)
    W1 = 0.421
    print("W1= ", W1)
    W2 = 0.927
    print("W2 = ", W2)
    kT = 0.501
    print("kT = ", kT)
    K_H_delta_B = fun_for_H_K_delta_B1(p, R, m1, m2, T1, T2, r0, W1, W2, kT, g)
    K1 = K_H_delta_B[0]
    H1 = K_H_delta_B[1]
    delta = K_H_delta_B[2]
    B1 = K_H_delta_B[3]
    W = P * (Cp - Cf) / (Cf - Cw)
    print("W = ", W)
    F = P * (Cp - Cw) / (Cf - Cw)
    print("F = ", F)
    omega0 = round(func_for_Ideal_Cascade_concentration_omega(Cf, Cp, Cw, H1, K1, P, W))
    print("omega0 = ", omega0)
    omega1 = round(0.75 * omega0)
    print("omega1 = ", omega1)
    Lambda_ideal = fun_for_Lambda_ideal(Cf, Cp, Cw, F, H1, K1, P, W)
    print("Lambda_ideal = ", Lambda_ideal)
    fun_for_B_opt_k_opt(Cf, Cp, Cw, F, H1, K1, P, W, delta, omega1)
    fun_for_plot_con_rectangle(K1, H1, Cw, Cp, Cf, P, W, omega1)
    return 0


if __name__ == '__main__':
    main()
