using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace Stirling
{
    internal class Program
    {
        static int nOParticles = 1000;
        static int step = 0;
        static double currentStepSize;
        static double globalTime = 0;

        static double[][] positions;
        static double[][] velocities;
        static double[] masses;
        static double particleSize = 40 * Math.Pow(10, 20/3 - 11); //0,05 was way too small - almost no particles collided

        static double smallCircleRadius = 0.06;
        static double circleRadius = 0.15;
        static double height = 1;

        static double hotPlate = 460 * Math.Pow(10, 21); //hotPlate is the upper one
        static double coldPlate = 200 * Math.Pow(10, 21);

        static double rOfFlyWheel = 0.1;
        static double upperPistonPos;
        static double lowerPistonPos;
        static double omega = 2;
        static double theta = 0.5 * 0.3 * 0.1 * 0.1;
        static double alpha = 0; //position of flywheel
        static double startingUpPistPosAngles = -Math.PI / 2; //0 angle is right - upper piston starts heading down

        static double k = 1.380649 * Math.Pow(10, -23);

        static int prevCollided1 = -10; //random number that cant be an index of anything - walls nor particles
        static int prevCollided2 = -10;


        static void Init()
        {
            Random rand = new Random();
            positions = new double[nOParticles][];
            velocities = new double[nOParticles][];
            masses = new double[nOParticles];

            for (int i = 0; i < nOParticles; i++)
            {
                positions[i] = Vector.Create(3);
                positions[i][0] = NextDouble(rand, 0 + particleSize, circleRadius - particleSize, 3);
                positions[i][1] = NextDouble(rand, 0 + particleSize, height - particleSize, 3);
                positions[i][2] = NextDouble(rand, 0, Math.Sqrt(Math.Pow(circleRadius, 2) - Math.Pow(positions[i][0], 2)), 3); //HELYRE KELL TENNI H NE CSAK A SZELERE MEHESSENEK

                velocities[i] = Vector.Create(3);
                velocities[i][0] = NextDouble(rand, -1, 1, 5);
                velocities[i][1] = NextDouble(rand, -1, 1, 5);
                velocities[i][2] = NextDouble(rand, -1, 1, 5);
                velocities[i] = Vector.MultiplyByNum(velocities[i], 60); //multiplying initial velocity like its in room temp

                masses[i] = 1.66 * Math.Pow(10, -5);
            }

            upperPistonPos = height + rOfFlyWheel * Math.Cos(alpha); //cos now but could be changed - with respect to the starting angledifference
            lowerPistonPos = rOfFlyWheel * Math.Sin(alpha);

        }

        static public double NextDouble(Random rand, double minValue, double maxValue, int decimalPlaces)
        {
            double randNumber = rand.NextDouble() * (maxValue - minValue) + minValue;
            return Convert.ToDouble(randNumber.ToString("f" + decimalPlaces));
        }

        static public double Particle_Particle_Collision(double[] v1, double[] v2, double[] p1, double[] p2, double r1, double r2)
        {
            double[] dpV = Vector.Substract(v1, v2);
            double[] p0V = Vector.Substract(p1, p2);
            double R = r1 + r2;

            double disP1 = Math.Pow(R, 2) * Vector.Square(dpV);
            double disP2 = Math.Pow(p0V[0] * dpV[1] - p0V[1] * dpV[0], 2);
            double disP3 = Math.Pow(p0V[1] * dpV[2] - p0V[2] * dpV[1], 2);
            double disP4 = Math.Pow(p0V[2] * dpV[0] - p0V[0] * dpV[2], 2);

            double firstPart = dpV[0] * p0V[0] + dpV[1] * p0V[1] + dpV[2] * p0V[2];
            double discriminant = disP1 - disP2 - disP3 - disP4;

            if (Vector.Square(dpV) == 0 || discriminant < 0)
            {
                return -1;
            }

            double t1 = (-firstPart + Math.Sqrt(discriminant)) / Vector.Square(dpV);
            double t2 = (-firstPart - Math.Sqrt(discriminant)) / Vector.Square(dpV);

            double t = (t1 > 0 && t2 > 0) ? Math.Min(t1, t2) : Math.Max(t1, t2);

            if (t < 0)
                return -1;
            else
                return t;
        }

        static public double Particle_Wall_Collision(double[] p, double[] v, double R, double r, double k, int wallType)
        {
            //walltype - -1 - cylinder palást - -2 - central circle with a hole - -3, -4 - pistons
            if (wallType == -1)
            {
                double radius = R - k;
                double temp = Math.Acos(p[0] / Vector.Magnitude(new double[] { p[0], p[2] }));
                double alpha = p[2] < 0 ? -temp : temp;
                temp = Math.Acos(v[0] / Vector.Magnitude(new double[] { v[0], v[2] }));
                double beta = v[2] < 0 ? -temp : temp;

                double multiplicant = Vector.Magnitude(new double[] { p[0], p[2] }) / Vector.Magnitude(new double[] { v[0], v[2] });
                double rootSecondPart = (radius / Vector.Magnitude(new double[] { p[0], p[2] }));
                double root = Math.Sqrt(Math.Pow(Math.Cos(alpha - beta), 2) + (rootSecondPart + 1) * (rootSecondPart - 1));
                double lastPart = Math.Cos(alpha - beta);


                double t1 = multiplicant * (root - lastPart);
                double t2 = multiplicant * -(root + lastPart);

                double t = (t1 > 0 && t2 > 0) ? Math.Min(t1, t2) : Math.Max(t1, t2);

                if (t < 0)
                    return -1;
                else
                    return t;

                //JOEZIGY?????
            }

            else if (wallType == -2) //kozepso talcaize
            {
                double t = (Math.Abs(height - p[1]) - k) / v[1];

                if (Math.Pow((p[0] + v[0] * t), 2) + Math.Pow((p[2] + v[2] * t), 2) < Math.Pow(r, 2) || t < 0)
                    return -1;
                else
                    return t;
            }

            else if (wallType == -3) //felso dugattyu
            {
                double t1 = (p[1] - k - upperPistonPos - rOfFlyWheel * Math.Cos(alpha)) / (-v[1] - rOfFlyWheel * Math.Sin(alpha) * omega);
                double t2 = (p[1] - k - upperPistonPos - rOfFlyWheel * Math.Cos(alpha)) / (-v[1] + rOfFlyWheel * Math.Sin(alpha) * omega);

                double t = (t1 > 0 && t2 > 0) ? Math.Min(t1, t2) : Math.Max(t1, t2);

                if (t < 0)
                    return -1;
                else
                    return t;
            }

            else if (wallType == -4) //also dugattyu
            {
                double t1 = (p[1] + k - lowerPistonPos - rOfFlyWheel * Math.Sin(alpha)) / (-v[1] - rOfFlyWheel * Math.Cos(alpha) * omega);
                double t2 = (p[1] + k - lowerPistonPos - rOfFlyWheel * Math.Sin(alpha)) / (-v[1] + rOfFlyWheel * Math.Cos(alpha) * omega);

                double t = (t1 > 0 && t2 > 0) ? Math.Min(t1, t2) : Math.Max(t1, t2);

                if (t < 0)
                    return -1;
                else
                    return t;
            }

            return -1;
        }

        static public void UpdatePositions(int n1, int n2)
        {
            alpha += omega * currentStepSize;

            upperPistonPos = height + rOfFlyWheel * Math.Cos(alpha);
            lowerPistonPos = rOfFlyWheel * Math.Sin(alpha);

            for (int i = 0; i < nOParticles; i++)
            {
                if (i == n1 && n2 >= 0)
                {
                    positions[n2][0] += velocities[n2][0] * currentStepSize;
                    positions[n2][1] += velocities[n2][1] * currentStepSize;
                    positions[n2][2] += velocities[n2][2] * currentStepSize;
                }
                else if (i == n1 && n2 < -2)
                {
                    positions[i][0] += velocities[i][0] * currentStepSize;
                    positions[i][1] = n2 == -3 ? upperPistonPos : lowerPistonPos;
                    positions[i][2] += velocities[i][2] * currentStepSize;
                    continue;
                }
                else if (i == n2)
                {
                    continue;
                }

                positions[i][0] += velocities[i][0] * currentStepSize;
                positions[i][1] += velocities[i][1] * currentStepSize;
                positions[i][2] += velocities[i][2] * currentStepSize;
            }
        }

        static public void UpdateVelocity(int nOParticle, int nOCollided)
        {
            if (nOCollided >= 0) //reszecskevel utkozott
            {
                double[] r0 = Vector.Substract(positions[(int)nOCollided], positions[nOParticle]);
                double[] r1 = Vector.Substract(positions[nOParticle], positions[(int)nOCollided]);
                double[] u0 = Vector.DivideByNum(r0, Vector.Magnitude(r0));
                double[] u1 = Vector.DivideByNum(r1, Vector.Magnitude(r1));
                double[] vi0 = Vector.MultiplyByNum(u0, (Vector.Magnitude(Vector.Add(r0, velocities[nOParticle])) / Vector.Magnitude(r0)));
                double[] vi1 = Vector.MultiplyByNum(u1, (Vector.Magnitude(Vector.Add(r1, velocities[(int)nOCollided])) / Vector.Magnitude(r1)));
                velocities[nOParticle] = Vector.Add(Vector.Substract(velocities[nOParticle], vi0), vi1);
                velocities[(int)nOCollided] = Vector.Add(Vector.Substract(velocities[(int)nOCollided], vi1), vi0);
            }
            else if (nOCollided == -1)//fallal utkozott
            {
                double[] pu = Vector.Add(positions[nOParticle], Vector.MultiplyByNum(velocities[nOParticle],
                    Particle_Wall_Collision(positions[nOParticle], velocities[nOParticle], circleRadius, smallCircleRadius, particleSize, -1)));
                double gamma = Math.PI / 2 - Math.Acos(Math.Abs(Vector.DotProduct(velocities[nOParticle], pu)) / (Vector.Magnitude(velocities[nOParticle]) * Vector.Magnitude(pu)));

                var temp = 3 * k * ((positions[nOParticle][1] >= height / 2) ? hotPlate : coldPlate) / masses[nOParticle];
                var d = Math.Sqrt(temp / Vector.Square(velocities[nOParticle]));

                velocities[nOParticle][0] = Math.Cos(2 * gamma) * velocities[nOParticle][0] - Math.Sin(2 * gamma) * velocities[nOParticle][2];
                velocities[nOParticle][2] = Math.Sin(2 * gamma) * velocities[nOParticle][0] + Math.Cos(2 * gamma) * velocities[nOParticle][2];
                velocities[nOParticle] = Vector.MultiplyByNum(velocities[nOParticle], d);
            }
            else if (nOCollided == -2)
            {
                velocities[nOParticle][1] = -velocities[nOParticle][1];
            }

            else if (nOCollided == -3) //felso dugattyu
            {
                double curSpeed = rOfFlyWheel * -Math.Sin(alpha);
                velocities[nOParticle][1] = 2 * curSpeed - velocities[nOParticle][1];

                omega = Math.Sqrt(2 * masses[nOParticle] * curSpeed * velocities[nOParticle][1] / theta + Math.Pow(omega, 2));
                //Console.WriteLine(velocities[nOParticle][0] + " " + velocities[nOParticle][1] + " " + velocities[nOParticle][2]);
            }

            else if (nOCollided == -4) //also dugattyu
            {
                double curSpeed = rOfFlyWheel * Math.Cos(alpha);
                velocities[nOParticle][1] = 2 * curSpeed - velocities[nOParticle][1];

                omega = Math.Sqrt(2 * masses[nOParticle] * curSpeed * velocities[nOParticle][1] / theta + Math.Pow(omega, 2));
                //Console.WriteLine(velocities[nOParticle][0] + " " + velocities[nOParticle][1] + " " + velocities[nOParticle][2]);
            }
        }

        static void Main(string[] args)
        {
            int[] forVel = new int[nOParticles * 150];

            for (int i = 0; i < forVel.Length; i++)
            {
                forVel[i] = i % 100;
            }

            Init();
            List<(double, int, int)> times = new List<(double, int, int)>();
            
                for (int i = 0; i < nOParticles; i++)
                {
                    if (i != nOParticles - 1)
                    {
                        for (int j = i + 1; j < nOParticles; j++)
                        {
                            times.Add((Math.Floor(Particle_Particle_Collision(velocities[i], velocities[j], positions[i], positions[j], particleSize, particleSize) * Math.Pow(10, 10)) / Math.Pow(10, 10), i, j));
                        }
                    }
                    for (int j = 0; j < 4; j++)
                    {
                        times.Add((Math.Floor(Particle_Wall_Collision(positions[i], velocities[i], circleRadius, smallCircleRadius, particleSize, -j - 1) * Math.Pow(10, 10)) / Math.Pow(10, 10), i, -j - 1));
                    }
                }

                int minI = 0;
                for (int i = 0; i < times.Count; i++)
                {
                    if (times[minI].Item1 == -1)
                    {
                        minI++;
                    }
                    if (times[i].Item1 < times[minI].Item1 && times[i].Item1 != -1 && times[i].Item1 > Math.Pow(10, -20)) //put in a hard cap in case of irrational numbers
                    {
                        minI = i;
                    }
                }
                currentStepSize = times[minI].Item1;
                step++;
                globalTime += currentStepSize;
                Console.WriteLine(times[minI].Item2 + " " + times[minI].Item3 + " " + currentStepSize);
                UpdatePositions(times[minI].Item2, times[minI].Item3);
                UpdateVelocity(times[minI].Item2, times[minI].Item3);
                prevCollided1 = times[minI].Item2;
                prevCollided2 = times[minI].Item3;

            while (true) //theres a problem!!!!!!!!!!
            {   
                //;rrowVBVOIbvoiwRV'OIwvoebv;oB;OUVBA;URVB;AURVBOd
                for (int i = times.Count - 1; i >= 0; i--)
                {
                    if (times[i].Item2 == prevCollided1 || times[i].Item3 == prevCollided1 || times[i].Item2 == prevCollided2 ||
                        (times[i].Item3 == prevCollided2 && (prevCollided2 >= 0 || prevCollided2 == -3 || prevCollided2 == -4)))
                        times.RemoveAt(i);
                    else
                        times[i] = (times[i].Item1 == -1 ? times[i].Item1 : times[i].Item1 - currentStepSize, times[i].Item2, times[i].Item3);
                }

                for (int i = 0; i < nOParticles; i++)
                {
                    if (i == prevCollided1 || i == prevCollided2)
                    {
                        continue;
                    }

                    times.Add((Math.Floor(Particle_Particle_Collision(velocities[prevCollided1], velocities[i], positions[prevCollided1], positions[i], particleSize, particleSize) *
                           Math.Pow(10, 10)) / Math.Pow(10, 10), prevCollided1, i));
                    
                    if (prevCollided2 >= 0)
                    {
                        times.Add((Math.Floor(Particle_Particle_Collision(velocities[prevCollided2], velocities[i], positions[prevCollided2], positions[i], particleSize, particleSize) *
                            Math.Pow(10, 10)) / Math.Pow(10, 10), prevCollided2, i));
                    }
                    else if (prevCollided2 == -3 || prevCollided2 == -4)
                    {
                        times.Add((Math.Floor(Particle_Wall_Collision(positions[i], velocities[i], circleRadius, smallCircleRadius, particleSize, prevCollided2) *
                            Math.Pow(10, 10)) / Math.Pow(10, 10), i, prevCollided2));
                    }
                }

                for (int i = -1; i >= -4; i--)
                {
                    times.Add((Math.Floor(Particle_Wall_Collision(positions[prevCollided1], velocities[prevCollided1], circleRadius, smallCircleRadius, particleSize, i) *
                        Math.Pow(10, 10)) / Math.Pow(10, 10), prevCollided1, i));

                    if (prevCollided2 >= 0)
                    {
                        times.Add((Math.Floor(Particle_Wall_Collision(positions[prevCollided2], velocities[prevCollided2], circleRadius, smallCircleRadius, particleSize, i) *
                            Math.Pow(10, 10)) / Math.Pow(10, 10), prevCollided2, i));
                    }
                }

                minI = 0;
                for (int i = 0; i < times.Count; i++)
                {
                    if (times[minI].Item1 == -1)
                    {
                        minI++;
                    }
                    if (times[i].Item1 < times[minI].Item1 && times[i].Item1 != -1 && times[i].Item1 > Math.Pow(10, -20)) //put in a hard cap in case of irrational numbers
                    {
                        minI = i;
                    }
                }
                currentStepSize = times[minI].Item1;
                step++;
                globalTime += currentStepSize;
                UpdatePositions(times[minI].Item2, times[minI].Item3);
                UpdateVelocity(times[minI].Item2, times[minI].Item3);
                prevCollided1 = times[minI].Item2;
                prevCollided2 = times[minI].Item3;
                //Console.WriteLine(times[minI].Item2 + " " + times[minI].Item3 + " " + currentStepSize);
                //Console.WriteLine(velocities[0][0] + " " + velocities[0][1] + " " + velocities[0][2]);

                ////////////////////////////////////////
                //Console.WriteLine(omega);

                double vatlfent = 0;
                double vatllent = 0;
                double vsqatlfent = 0;
                double vsqatllent = 0;
                double noFent = 0;
                double noLent = 0;

                for (int i = 0; i < nOParticles; i++)
                {
                    if (positions[i][1] >= height / 2)
                    {
                        vatlfent += Vector.Magnitude(velocities[i]);
                        vsqatlfent += Vector.Square(velocities[i]);
                        noFent++;
                    }
                    else
                    {
                        vatllent += Vector.Magnitude(velocities[i]);
                        vsqatllent += Vector.Square(velocities[i]);
                        noLent++;
                    }
                }

                vatlfent /= nOParticles;
                vatllent /= nOParticles;
                vsqatlfent /= nOParticles;
                vsqatllent /= nOParticles;

                double Vfelso = (upperPistonPos - height/2) * Math.Pow(circleRadius, 2) * Math.PI;
                double Valso = (height/2 - lowerPistonPos) * Math.Pow(circleRadius, 2) * Math.PI;

                double rofent = masses[0] * noFent / Vfelso;
                double rolent = masses[0] * noLent / Valso;

                using (StreamWriter writetext = new StreamWriter("e:\\Developing\\TDK-stirling\\omega.txt", true))
                {
                    writetext.WriteLine(omega);
                }

                using (StreamWriter writetext = new StreamWriter("e:\\Developing\\TDK-stirling\\alpha.txt", true))
                {
                    writetext.WriteLine(alpha);
                }

                using (StreamWriter writetext = new StreamWriter("e:\\Developing\\TDK-stirling\\dt-global.txt", true))
                {
                    writetext.WriteLine(globalTime);
                }

                using (StreamWriter writetext = new StreamWriter("e:\\Developing\\TDK-stirling\\pfelso.txt", true))
                {
                    writetext.WriteLine(rofent * vatlfent);
                }

                using (StreamWriter writetext = new StreamWriter("e:\\Developing\\TDK-stirling\\palso.txt", true))
                {
                    writetext.WriteLine(rolent * vatllent);
                }

                using (StreamWriter writetext = new StreamWriter("e:\\Developing\\TDK-stirling\\Tfelso.txt", true))
                {
                    writetext.WriteLine(vsqatlfent);
                }

                using (StreamWriter writetext = new StreamWriter("e:\\Developing\\TDK-stirling\\Talso.txt", true))
                {
                    writetext.WriteLine(vsqatllent);
                }

                using (StreamWriter writetext = new StreamWriter("e:\\Developing\\TDK-stirling\\Vfelso.txt", true))
                {
                    writetext.WriteLine(Vfelso);
                }

                using (StreamWriter writetext = new StreamWriter("e:\\Developing\\TDK-stirling\\Valso.txt", true))
                {
                    writetext.WriteLine(Valso);
                }

                if (forVel.Contains(prevCollided1))
                {
                    forVel[Array.IndexOf(forVel, prevCollided1)] = -1;
                    if (prevCollided2 > 0 && forVel.Contains(prevCollided2))
                    {
                        forVel[Array.IndexOf(forVel, prevCollided2)] = -1;
                    }
                }

                Console.WriteLine(globalTime);

                /*if (forVel.All(x => x == -1))
                {
                    for (int i = 0; i < velocities.GetLength(0); i++)
                    {
                        Console.WriteLine(Vector.Magnitude(velocities[i]));
                    }
                    break;
                }*/

                ////////////////////////////////////////////
            }
        }
    }

    internal class Vector
    {
        public static double[] Create(int rows)
        {
            return new double[rows];
        }

        public static double[] CreateZero(int rows)
        {
            var temp = new double[rows];
            for (int i = 0; i < temp.Length; i++)
            {
                temp[i] = 0;
            }
            return temp;
        }

        public static double[] Substract(double[] v1, double[] v2)
        {
            if (v1.Length == v2.Length)
            {
                int l = v1.Length;
                double[] temp = new double[l];
                for (int i = 0; i < l; i++)
                {
                    temp[i] = v1[i] - v2[i];
                }
                return temp;
            }
            return new double[] { -1 };
        }

        public static double[] Add(double[] v1, double[] v2)
        {
            if (v1.Length == v2.Length)
            {
                int l = v1.Length;
                double[] temp = new double[l];
                for (int i = 0; i < l; i++)
                {
                    temp[i] = v1[i] + v2[i];
                }
                return temp;
            }
            return new double[] { -1 };
        }

        public static double[] DivideByNum(double[] v, double num)
        {
            int l = v.Length;
            double[] temp = new double[l];
            for (int i = 0; i < l; i++)
            {
                temp[i] = v[i] / num;
            }
            return temp;
        }

        public static double[] MultiplyByNum(double[] v, double num)
        {
            int l = v.Length;
            double[] temp = new double[l];
            for (int i = 0; i < l; i++)
            {
                temp[i] = v[i] * num;
            }
            return temp;
        }

        public static double DotProduct(double[] v1, double[] v2)
        {
            double temp = 0;
            for (int i = 0; i < v1.Length; i++)
            {
                temp += v1[i] * v2[i];
            }
            return Math.Abs(temp);
        }

        public static double DotProduct(double v1L, double v2L, double theta)
        {
            return v1L * v2L * theta;
        }

        public static double Magnitude(double[] v)
        {
            if (v.Length == 0)
            {
                return -1;
            }
            double temp = 0;
            for (int i = 0; i < v.Length; i++)
            {
                temp += v[i] * v[i];
            }
            return Math.Sqrt(temp);
        }

        public static double Square(double[] v)
        {
            return Math.Pow(Magnitude(v), 2);
        }

        public static double[] UnitV(double[] v)
        {
            double[] temp = new double[v.Length];
            for (int i = 0; i < v.Length; i++)
            {
                temp[i] = v[i] / Vector.Magnitude(v);
            }

            return temp;
        }
    }
}
