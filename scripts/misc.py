from matplotlib.pyplot import plot as plt

def booth(solution: tuple) -> float:
    """Booth function

    Args:
        solution (tuple): (x, y) solution to the function

    Returns:
        float: z value of solution
    """
    x,y = solution
    
    z = pow((x + 2*y - 7), 2) + pow((2*x + y - 5), 2)

def bohachevsky_one(solution: tuple) -> float:
    """Bohachevsky1 objective function.

    Args:
        solution (tuple): (x, y) solution to the function

    Returns:
        float: z value of solution
    """
if __name__ == "__main__":
    for i in range(1,17):
        v = pow(3,i) - pow(2,i) - 1
        print("solution length: ", i, " | num schemata: ", v)