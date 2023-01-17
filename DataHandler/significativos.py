"""
Author: Thais Silva Abelha
Email: thais.silva.abelha@gmail.com

Coloque o valor e a incerteza (svalor), para formatá-los com 2 algarismos
significativos
"""


class dois_significativos:
    def __init__(self, valor: float, svalor: float) -> None:
        self.valor = valor
        self.svalor = svalor

    def incerteza(self) -> tuple:
        s = "%.2g" % (self.svalor)
        i = s.find(".")
        f = len(s)

        if i != -1:
            n = f - i - 1
            if n == 1 and self.svalor < 1:
                n = n + 1
            decimais = "%0." + str(n) + "f"
            svalor_final = decimais % (self.svalor)
        if i == -1:
            if len(s) == 1:
                n = 1
                decimais = "%0." + str(n) + "f"
                svalor_final = decimais % (self.svalor)
            else:
                n = 0
                decimais = "%0." + str(n) + "f"
                svalor_final = decimais % (self.svalor)

        return svalor_final, decimais, n

    def forma1(self) -> str:
        svalor_final, decimais, n = self.incerteza()
        return decimais % self.valor + "#" + svalor_final

    def forma2(self) -> str:
        svalor_final, decimais, n = self.incerteza()
        return (
            decimais % self.valor
            + "("
            + "%1.0f" % (float(svalor_final) * 10**n)
            + ")"
        )
