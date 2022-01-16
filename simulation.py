# Klassenmodellierung des SEIR-Modell
# Floyd Erdmann

# Import des Moduls 'Numpy' zur Berechnung
import numpy as np

# Import des Moduls 'Matplotlib' zur Darstellung
import matplotlib.pyplot as plt

class SEIR:
    """
        - die initialen Werte entsprechen den Kategorien S, E, I und R
        - die dynamischen Paramter alpha, beta, gamma und rho sind im 
          internen Array 'parameter_' abgelegt
    """
    def __init__(
        self,
        initiale_werte = [1 - 1/1000, 1/1000, 0, 0],
        parameter_     = [0.2, 1.75, 0.5, 0.1]):

        # initiale Werte
        self.s0 = initiale_werte[0]
        self.e0 = initiale_werte[1]
        self.i0 = initiale_werte[2]
        self.r0 = initiale_werte[3]

        # Listen
        self.s = [self.s0]
        self.e = [self.e0]
        self.i = [self.i0]
        self.r = [self.r0]

        # dynamische Parameter
        self.alpha = parameter_[0]  # Inkubationszeit
        self.beta  = parameter_[0]  # Infektionsrate
        self.gamma = parameter_[0]  # Erholungszeit
        self.rho   = parameter_[0]  # Social distancing

        # Alle Parameter in einer Liste
        self.parameter_ = [self.alpha, self.beta, self.gamma, self.rho]

        # Alle finalen Werte in einer Liste
        self.werte_ = [self.s[-1], self.e[-1], self.i[-1], self.r[-1]]

    """
        - erneute Initialisierung mit neuen Werten
        - Parameter verbose: wenn 'True' werden die Werte ausgegeben
    """
    def reinitialisierung(self, initiale_werte, verbose=True):
        # Absicherung der Anzahl der Werte und des Datentyps
        assert len(initiale_werte)  ==4,    "Vier initiale Werte sind notwendig"
        assert type(initiale_werte) ==list, "Initiale Werte werden als Liste erwartet"

        # Initiale Werte
        self.s0 = initiale_werte[0]
        self.e0 = initiale_werte[1]
        self.i0 = initiale_werte[2]
        self.r0 = initiale_werte[3]

        if verbose:
            print("Initialisierung erfolgt mit den Werten:\n")
            print("S0: ",self.s0)
            print("E0: ",self.e0)
            print("I0: ",self.i0)
            print("R0: ",self.r0)

    """
        - setzen der dynamischen Parameter
    """
    def setze_parameter(self, parameter_, verbose=True):
        # Absicherung der Anzahl der Parameter und des Datentyps
        assert len(parameter_)  == 4,    "Four parameter values are expected"
        assert type(parameter_) == list, "Parameter values are expected in a list"
        
        # dynamische Parameter
        self.alpha      = parameter_[0]
        self.beta       = parameter_[1]
        self.gamma      = parameter_[2]
        self.rho        = parameter_[3]
        self.parameter_ = [self.alpha,self.beta,self.gamma,self.rho]
        
        if verbose:
            print("Setzen der Parameter mit den folgenden Werten:\n")
            print("alpha: ",self.alpha)
            print("beta: " ,self.beta)
            print("gamma: ",self.gamma)
            print("rho: "  ,self.rho)

    """
        - zurücksetzen der internen Liste in den Ausgangsstauts ('zero-state')
    """
    def zuruecksetzen(self):
        self.s = [self.s0]
        self.e = [self.e0]
        self.i = [self.i0]
        self.r = [self.r0]

    """
        - starten der dynamischen Simulation
        - Argumente:
                - zeit_max: maximale Simulationsdauer, z.B. 20 oder 100 (als Tage)
                - zeit_intervall: Zeitintervall, z.B. 0.1 or 0.02
                - zuruecksetzen: A flag to zuruecksetzen the internal lists (restarts the simulation from initial values)
    """
    def start(self, zeit_max=100, zeit_intervall=0.1, zuruecksetzen=True):
        
        if zuruecksetzen:
            self.zuruecksetzen()
        
        # Array für Zeitintervalle
        t = np.linspace(0, zeit_max, int(zeit_max/zeit_intervall) + 1)
        
        # temporäre Listen
        S = self.s
        E = self.e
        I = self.i
        R = self.r
        
        # temporäre Parameter
        alpha  = self.alpha
        beta   = self.beta
        gamma  = self.gamma
        rho    = self.rho
        zeit_intervall = t[1] - t[0]
        
        # Schleife zur Berechnung (siehe SEIR-Modell)
        for _ in t[1:]:
            next_S = S[-1] - ( rho * beta * S[-1] * I[-1] ) * zeit_intervall
            next_E = E[-1] + ( rho * beta * S[-1] * I[-1] - alpha * E[-1]) * zeit_intervall
            next_I = I[-1] + ( alpha * E[-1] - gamma * I[-1] ) * zeit_intervall
            next_R = R[-1] + ( gamma * I[-1]) * zeit_intervall
            S.append(next_S)
            E.append(next_E)
            I.append(next_I)
            R.append(next_R)
        
        # 'Sammeln' der Ergebniss
        ergebnis = np.stack([S, E, I, R]).T
        self.s, self.e, self.i, self.r = S, E, I, R
        
        # Schreiben der endgültigen Werte 
        self.werte_ = [self.s[-1], self.e[-1], self.i[-1], self.r[-1]]
        
        return ergebnis


################# Zeichenfunktionen #########################################

    
    """
        Zeichnen der Ergebnisse einfachen Ergebnisse
    """
    def zeichnen(self, ergebnis=None):
        # Startet eine Simulation wenn kein Ergebniss aufgerufen wird
        if ergebnis is None:
            ergebnis = self.run()
        
        # Zeichnen
        plt.figure(figsize = (12,8))
        plt.plot(ergebnis, lw=3)
        plt.title('SEIR-Modell', fontsize=18)
        plt.legend(['Anfällige', 'Ausgesetzte', 'Infizierte', 'Genesene'], fontsize=15)
        plt.xlabel('Zeitintervall', fontsize=16)
        plt.ylabel('Anteil der Bevölkerung', fontsize=16)
        plt.grid(True)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.show()


    """
        Plots the given variable
        Expect a list or Numpy array as the variable
        If var is None, plots the infected fraction
    """
    def zeichnen_variable( self, var, var_name=None, show=True):
        if var is None:
            var = self.i
        plt.figure(figsize=(12,8))
        plt.plot(var,lw=3,c='blue')
        plt.title('Basic SEIR Model',fontsize=18)
        if var_name is not None:
            plt.legend([var_name],fontsize=15)
        plt.xlabel('Time Steps',fontsize=16)
        plt.ylabel('Fraction of Population',fontsize=16)
        plt.grid(True)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        if show:
            plt.show()











