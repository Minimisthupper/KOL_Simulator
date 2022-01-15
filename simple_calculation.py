from numpy import array as vector

# Explizites Euler-Verfahren
def euler_method(f, t0, x0, t1, h):
    t = t0; x = x0
    a = [[t, x]]
    for k in range(0, 1 + int((t1 - t0)/h)):
        t = t0 + k*h
        x = x + h*f(t, x)
        a.append([t, x])
    return a

def SEIR_model(alpha, beta, gamma):
    def f(t, x):
        s, e, i, r = x
        return vector([
            -beta*s*i,
            beta*s*i - alpha*e,
            alpha*e - gamma*i,
            gamma*i
        ])
    return f

def SEIR_simulation(alpha, beta, gamma, e0, i0, days, step=0.1):
    x0 = vector([1.0 - e0 - i0, e0, i0, 0.0])
    f = SEIR_model(alpha, beta, gamma)
    return euler_method(f, 0, x0, days, step)

def diagram(simulation):
    import matplotlib.pyplot as plot
    # plot.style.use('fivethirtyeight')
    figure,axes = plot.subplots()
    figure.subplots_adjust(bottom = 0.15)
    axes.grid(linestyle = ':', linewidth = 2.0, color = "#808080")
    t,x = zip(*simulation())
    s, e, i, r = zip(*x)
    axes.plot(t, s, color = "#0000cc")
    axes.plot(t, e, color = "#ffb000", linestyle = '--')
    axes.plot(t, i, color = "#a00060")
    axes.plot(t, r, color = "#008000", linestyle = '--')
    plot.show()

def simulation1():
    N = 83200000 # Einwohnerzahl von Deutschland 2019/2020
    R0 = 2.4; gamma = 1/3.0
    return SEIR_simulation(
        alpha = 1/5.5, beta = R0*gamma, gamma = gamma,
        e0 = 40000.0/N, i0 = 10000.0/N, days = 140)

diagram(simulation1)