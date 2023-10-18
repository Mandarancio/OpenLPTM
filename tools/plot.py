import sys

import matplotlib.pyplot as plt
import pandas as pd

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("missing CSV file path")
        sys.exit(1)
    x = pd.read_csv(sys.argv[1], header=None, comment="#").values
    plt.plot(x)
    plt.show()
