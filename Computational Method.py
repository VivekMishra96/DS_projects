class Root_Finding:

    class Newton_Raphson:
        def func(x):
            return x * x * x - x * x + 2

        def derivFunc(x):
            return 3 * x * x - 2 * x

        def newtonRaphson():
            x = int(input("Enter your guess value: "))
            h = Root_Finding.Newton_Raphson.func(
                x) / Root_Finding.Newton_Raphson.derivFunc(x)
            while abs(h) >= 0.0001:
                h = Root_Finding.Newton_Raphson.func(
                    x) / Root_Finding.Newton_Raphson.derivFunc(x)
                x = x - h
            print("The value of the root is : ", "%.4f" % x)

    class Fixed_Point_Iteration:
        def f(x):
            return x*x*x + x*x - 1

        def g(x):
            import math
            return 1/math.sqrt(1+x)

        def fixed_point_iteration():
            x0 = float(input('Enter Guess: '))
            e = float(input('Tolerable Error: '))
            N = int(input('Maximum Step: '))
            step = 1
            flag = 1
            condition = True
            while condition:
                x1 = Root_Finding.Fixed_Point_Iteration.g(x0)
                print('Iteration-%d, x1 = %0.6f and f(x1) = %0.6f' %
                      (step, x1, Root_Finding.Fixed_Point_Iteration.f(x1)))
                x0 = x1
                step = step + 1

                if step > N:
                    flag = 0
                    break

                condition = abs(Root_Finding.Fixed_Point_Iteration.f(x1)) > e

            if flag == 1:
                print('\nRequired root is: %0.8f' % x1)
            else:
                print('\nNot Convergent.')

    class Secant_Method:
        def f(x):
            f = pow(x, 3) + x - 1
            return f

        def secant_method():
            x1 = int(input("Enter lower limit: "))
            x2 = int(input("Enter upper limit: "))
            E = float(input("Enter tolerence: "))
            n = 0
            xm = 0
            x0 = 0
            c = 0
            if (Root_Finding.Secant_Method.f(x1) * Root_Finding.Secant_Method.f(x2) < 0):
                while True:

                    x0 = ((x1 * Root_Finding.Secant_Method.f(x2) - x2 * Root_Finding.Secant_Method.f(x1)) /
                          (Root_Finding.Secant_Method.f(x2) - Root_Finding.Secant_Method.f(x1)))
                    c = Root_Finding.Secant_Method.f(
                        x1) * Root_Finding.Secant_Method.f(x0)
                    x1 = x2
                    x2 = x0
                    n += 1

                    if (c == 0):
                        break
                    xm = ((x1 * Root_Finding.Secant_Method.f(x2) - x2 * Root_Finding.Secant_Method.f(x1)) /
                          (Root_Finding.Secant_Method.f(x2) - Root_Finding.Secant_Method.f(x1)))

                    if(abs(xm - x0) < E):
                        break

                print("Root of the given equation =",
                      round(x0, 6))
                print("No. of iterations = ", n)

            else:
                print("Can not find a root in ",
                      "the given interval")

    class False_Position_Method:
        def f(x):
            return x**3-5*x-9

        def false_position_method():
            x0 = int(input("Enter lower limit: "))
            x1 = int(input("Enter upper limit: "))
            e = float(input("Enter tolerence: "))
            step = 1
            condition = True
            while condition:
                x2 = x0 - (x1-x0) * Root_Finding.False_Position_Method.f(x0)/(
                    Root_Finding.False_Position_Method.f(x1) - Root_Finding.False_Position_Method.f(x0))
                print('Iteration-%d, x2 = %0.6f and f(x2) = %0.6f' %
                      (step, x2, Root_Finding.False_Position_Method.f(x2)))

                if Root_Finding.False_Position_Method.f(x0) * Root_Finding.False_Position_Method.f(x2) < 0:
                    x1 = x2
                else:
                    x0 = x2

                step = step + 1
                condition = abs(Root_Finding.False_Position_Method.f(x2)) > e

            print('\nRequired root is: %0.8f' % x2)


class Gauss_Elimination:
    def gauss_elimination():
        import numpy as np
        import sys

        n = int(input('Enter number of unknowns: '))
        a = np.zeros((n, n+1))
        x = np.zeros(n)
        print('Enter Augmented Matrix Coefficients:')
        for i in range(n):
            for j in range(n+1):
                a[i][j] = float(input('a['+str(i)+'][' + str(j)+']='))

        for i in range(n):
            if a[i][i] == 0.0:
                sys.exit('Divide by zero detected!')

            for j in range(i+1, n):
                ratio = a[j][i]/a[i][i]

                for k in range(n+1):
                    a[j][k] = a[j][k] - ratio * a[i][k]

        x[n-1] = a[n-1][n]/a[n-1][n-1]

        for i in range(n-2, -1, -1):
            x[i] = a[i][n]

            for j in range(i+1, n):
                x[i] = x[i] - a[i][j]*x[j]

            x[i] = x[i]/a[i][i]

        print('\nRequired solution is: ')
        for i in range(n):
            print('X%d = %0.2f' % (i, x[i]), end='\t')


class Lagrange_Interpolation:
    def lagrange_interpolation():
        import numpy as np
        n = int(input('Enter number of data points: '))
        x = np.zeros((n))
        y = np.zeros((n))
        print('Enter data for x and y: ')
        for i in range(n):
            x[i] = float(input('x['+str(i)+']='))
            y[i] = float(input('y['+str(i)+']='))

        xp = float(input('Enter interpolation point: '))
        yp = 0
        for i in range(n):

            p = 1

            for j in range(n):
                if i != j:
                    p = p * (xp - x[j])/(x[i] - x[j])

            yp = yp + p * y[i]

        print('Interpolated value at %.3f is %.3f.' % (xp, yp))


class Numerical_Integration:
    class Trapezoidal_Method:
        def f(x):
            return 1/(1 + x**2)

        def trapezoidal_method():
            x0 = float(input("Enter lower limit: "))
            xn = float(input("Enter upper limit: "))
            n = int(input("Enter no of sub intervals: "))
            h = (xn - x0) / n
            integration = Numerical_Integration.Trapezoidal_Method.f(
                x0) + Numerical_Integration.Trapezoidal_Method.f(xn)

            for i in range(1, n):
                k = x0 + i*h
                integration = integration + 2 * \
                    Numerical_Integration.Trapezoidal_Method.f(k)

            integration = integration * h/2
            print("Integration result by Trapezoidal method is: %0.6f" %
                  (integration))

    class Simpson_1_3_Method:
        def f(x):
            return 1/(1 + x**2)

        def simpson_method():
            x0 = float(input("Enter lower limit: "))
            xn = float(input("Enter higher limit: "))
            n = int(input("Enter no of sub intervals: "))
            h = (xn - x0) / n

            integration = Numerical_Integration.Simpson_1_3_Method.f(
                x0) + Numerical_Integration.Simpson_1_3_Method.f(xn)

            for i in range(1, n):
                k = x0 + i*h

                if i % 2 == 0:
                    integration = integration + 2 * \
                        Numerical_Integration.Simpson_1_3_Method.f(k)
                else:
                    integration = integration + 4 * \
                        Numerical_Integration.Simpson_1_3_Method.f(k)

            integration = integration * h/3
            print("Integration result by Trapezoidal method is: %0.6f" %
                  (integration))

    class Simpson_3_8_Method:
        def f(x):
            return 1/(1 + x**2)

        def simpson_3_8():
            x0 = float(input("Enter lower limit: "))
            xn = float(input("Enter higher limit: "))
            n = int(input("Enter no of sub intervals: "))
            h = (xn - x0) / n
            integration = Numerical_Integration.Simpson_3_8_Method.f(
                x0) + Numerical_Integration.Simpson_3_8_Method.f(xn)

            for i in range(1, n):
                k = x0 + i*h

                if i % 2 == 0:
                    integration = integration + 2 * \
                        Numerical_Integration.Simpson_3_8_Method.f(k)
                else:
                    integration = integration + 3 * \
                        Numerical_Integration.Simpson_3_8_Method.f(k)

            integration = integration * 3 * h / 8

            print("Integration result by Trapezoidal method is: %0.6f" %
                  (integration))


class ODE:
    class Euler_Method:
        def f(x, y):
            return (x + y + x * y)

        def euler_method():
            x0 = float(input("Enter lower limit: "))
            y0 = float(input("Enter higher limit: "))
            x = float(input("Enter the value at which you need approximation: "))
            n = int(input("Enter no of sub intervals: "))
            h = (x - x0) / n
            temp = -0
            while x0 < x:
                temp = y0
                y0 = y0 + h * ODE.Euler_Method.f(x0, y0)
                x0 = x0 + h

            print("Approximate solution at x = ", x, " is ", "%.6f" % y0)
ODE.Euler_Method.euler_method()