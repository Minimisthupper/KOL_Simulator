# KOL, Floyd Erdmann
# Klassenmodellierung des SEIR-Modell
# 
# Fragestellung:
#   - Wie wirkt sich die Globalisierung auf die Verbreitung von Krankheitserregern aus.

# Import des Moduls 'Numpy' zur Berechnung
import numpy as np

# Import des Moduls 'Matplotlib' zur Darstellung
import matplotlib.pyplot as plt


"""
    - die initialen Werte entsprechen den Kategorien S, E, I und R
    - die dynamischen Paramter alpha, beta, gamma und rho sind im 
      internen Array 'parameter_' abgelegt
    Funktions-Argumente:
        - initiale_werte:
            - Anteil der Bevoelkerung in den Kategorien S, E, I und R
            - Zu Beginn der Ausbreitung sind nur 0.001% Bevoelkerung Ausgesetzt.
"""
class SEIR:
    """
        - __init__ ist eine 'BuildIn'-Funktion
        - wird ausgeführt um die Klasse zu initialisieren (ein Objekt zu erstellen)
        - dient der Zuweisung don Werten zu Objekt-Attributen

        - Initiale parameters (in 'params_' gespeichert)
            alpha (1 / Inkubationszeit)                 = 0.2
            beta  (durchschnittliche Kontaktrate)       = 1.75
            gamma ( 1 / mittlere Zeit der Regeneration) = 0.5
            rho   (Social Distancing Faktor)            = 0.1

        - rho = 1 bedeutet "komplettes Social Distancing")

        - Quelle für reale Covid-Werte aus dem Jahre 2019:
        https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(20)30074-7/fulltext
    """
    def __init__(
        self,
        # Anfällige, Ausgesetzte, Infizierte, Genesene
        initiale_werte = [1 - 1/1000, 1/1000, 0, 0],
        # Inkubatioszeit, Kontaktrate, Regeneration, Social Distancing
        parameter_     = [0.2, 1.80, 0.5, 1.0]):
        

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
        self.alpha = parameter_[0]
        self.beta  = parameter_[1]
        self.gamma = parameter_[2]
        self.rho   = parameter_[3]
        
        # Alle Parameter in einer Liste
        self.parameter_ = [self.alpha,self.beta,self.gamma,self.rho]
        
        # Alle finalen Werte in einer Liste
        self.vals_ = [self.s[-1], self.e[-1], self.i[-1], self.r[-1]]
    
    """
        - erneute Initialisierung mit neuen Werten
        - Parameter 'verbose': wenn 'True' werden die Werte ausgegeben
    """
    def reinitialisierung(
        self,
        initiale_werte,
        verbose=False):
    
        assert len(initiale_werte)  == 4,    "Vier initiale Werte sind notwendig"
        assert type(initiale_werte) == list, "Initiale Werte werden als Liste erwartet"
        
        # Initial values
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
    def setze_parameter(
        self,
        parameter_,
        verbose=False):
        
        assert len(parameter_)  == 4,    "Vier Parameter-Werte sind notwendig"
        assert type(parameter_) == list, "Parameter-Werte werden als Liste erwartet"
        
        # dynamische Parameter
        self.alpha      = parameter_[0]
        self.beta       = parameter_[1]
        self.gamma      = parameter_[2]
        self.rho        = parameter_[3]
        self.parameter_ = [self.alpha,self.beta,self.gamma,self.rho]
        
        if verbose:
            print("Setzen der Parameter mit den folgenden Werten:\n")
            print("alpha: ",self.alpha)
            print("beta:  ",self.beta)
            print("gamma: ",self.gamma)
            print("rho:   ",self.rho)

    
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
    def start(
        self,
        zeit_max         = 100,
        zeit_intervall   = 0.1,
        zuruecksetzen    = True):
        
        if zuruecksetzen:
            self.zuruecksetzen()
        
        
        """
            Dokumentation zur Funktion 'linspace':
                - https://numpy.org/doc/stable/reference/generated/numpy.linspace.html
                - Gibt Zahlen in gleichmäßigen Abstand über ein spezifiziertes Intervall aus.
        """
        # Array für Zeitintervalle
        t = np.linspace(0, zeit_max, int( zeit_max / zeit_intervall )  + 1)
        
        # temporäre Listen
        S = self.s
        E = self.e
        I = self.i
        R = self.r
        
        # temporäre Parameter
        alpha = self.alpha
        beta  = self.beta
        gamma = self.gamma
        rho   = self.rho
        zeit_intervall    = t[1] - t[0]
        
        # Schleife zur Berechnung (siehe Kapitel SEIR-Modell)
        for _ in t[1:]:
            next_S = S[-1] - ( rho * beta *S[-1] * I[-1] ) * zeit_intervall
            next_E = E[-1] + ( rho * beta* S[-1] * I[-1] - alpha * E[-1]) * zeit_intervall
            next_I = I[-1] + ( alpha * E[-1] - gamma*I[-1]) * zeit_intervall
            next_R = R[-1] + ( gamma * I[-1]) * zeit_intervall
            
            S.append(next_S)
            E.append(next_E)
            I.append(next_I)
            R.append(next_R)
        
        """
            Dokumentation zur Funktion 'linspace':
                - https://numpy.org/doc/stable/reference/generated/numpy.stack
                - Fügt die Arrays zusammen.
        """
        # 'Einsammeln' der Ergebnisse
        ergebnis = np.stack([S, E, I, R]).T
        self.s = S
        self.e = E
        self.i = I
        self.r = R
        
        # Schreiben der finalen Werte in das Array
        self.vals_ = [self.s[-1], self.e[-1], self.i[-1], self.r[-1]]
        
        return ergebnis


################# Zeichenfunktionen #########################################
    
    """
        Zeichnen der Ergebnisse mit Hilfe des Moduls Matplotlib

        Dokumentation: https://matplotlib.org/stable/api/index.html
    """
    def zeichnen(
        self,
        ergebnisse = None):
        
        # Runs a simulation is no result is provided
        if ergebnisse is None:
            ergebnisse = self.start()
        
        # Darstellung
        # Top-Level Widget: beinhaltet alle Zeichenelemente
        """
        class matplotlib.figure.Figure(
            figsize=None, dpi=None, facecolor=None, edgecolor=None, linewidth=0.0,
            frameon=None, subplotpars=None, tight_layout=None, constrained_layout=None, *, layout=None, **kwargs)
        """
        plt.figure(figsize=(15,11), linewidth=3)
        
        # matplotlib.pyplot.plot(*args, scalex=True, scaley=True, data=None, **kwargs)
        plt.plot(ergebnisse, linewidth=3)

        # matplotlib.pyplot.title(label, fontdict=None, loc=None, pad=None, *, y=None, **kwargs)
        plt.title(
            label='Simulation der Ausbreitung von Krankheitserregern\n (erweitertes SEIR-Modell)',
            loc='center',
            fontsize=28)
        
        # Legende
        plt.legend(
            ['Anfällige (Susceptible)', 'Ausgesetzte (Exposed)', 'Infizierte (Infected)', 'Genesene (Recovered)'],
            fontsize=20)

        # x-Achse
        plt.xlabel(
            'Zeitintervall (in Tagen)',
            fontsize=20)

        # y-Achse
        plt.ylabel(
            'Anteil der Bevölkerung (in %)',
            fontsize=20)
        
        # Gridanzeige
        #plt.grid(True)
        plt.grid(linestyle='dashdot', linewidth=0.5)

        # Bezeichner/Label der x- und y-Achse
        plt.xticks(fontsize=15);
        plt.yticks(fontsize=15)
        
        # Diagramm anzeigen
        plt.show()