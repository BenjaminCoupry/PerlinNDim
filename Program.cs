using System;
using System.Drawing;
using System.Windows.Media.Imaging;
using System.Threading;
using System.Collections.Generic;

namespace Geologik
{
    class Program
    {
        public enum TypeInterpolation { Lineaire, Cosinus, Hermite, C2 };
        static void Main(string[] args)
        {
            Random r = new Random(454542);
            double scale = 0.005;
            Application p1 = perlin(10000, TypeInterpolation.C2, 7, 0.1, 0.5, 5, 1.1, 0, 0.7);
            Application shift = (x) => ((p1(x) * 2) - 1);
            Application p2 = perlin(55656, TypeInterpolation.C2, 7, 3, 0.4, 0.5, 1.1, 0, 0.7);
            double[,] crateres = sample2D(300, 300, 1.0 / 100, p1);
            crateres = Bombarder(crateres,20, 1,ref r, shift,15, 110, 0.1, 1, 0.1, 1, 0.05, 1,5,0.5);
            visu2D(crateres, false).Save("D:/lab/Terrain/Crateres.bmp");
            /**
            Application rings = Ondelettes(2, new double[] { -0.5, -0.5 });
            double[,] grille = sample2D(100, 100, 1.0/100, p1);
            Normer(grille);
            visu2D(grille, false).Save("D:/lab/Terrain/_perlinPur.bmp");
            //grille = froisser(grille,new double[] {0.5,0.5 },new double[] { 0,0},0.3,2,Cassure(200),p2,p1,scale,1,0.3,1);
            //visu2D(grille, false).Save("D:/lab/Terrain/_perlinGeol.bmp");
            
            for (int i = 0; i < 5; i++)
             {
                Console.WriteLine(i);
                erosionTerrainCroisee(ref r, ref grille, new double[] { 0.2, 1 }, 600, 0.5, 0.0002, 3, 0.05, 0.05, 10, scale);
                grille = lissage(grille, 1);
                visu2D(grille, false).Save("D:/lab/Terrain/historique/perlinErosion"+i+".bmp");
                exportCsv(grille, "D:/lab/Terrain/historique/perlinErosion" + i + ".csv");
            }
            double[,] lacs = altitudesLacs(grille, 0, 0.005);
            visu2D(lacs, false).Save("D:/lab/Terrain/LACS.bmp");
            **/
            /**
            Frame fr = () => { erosionTerrainCroisee(ref r, ref grille, new double[] { 0.2, 1 }, 1000, 0.5, 0.0002, 3, 0.05, 0.05, 2600, scale);
                return visu2D(grille, false); };
            creerGIF("D:/lab/Terrain/historique/perlinErosion.gif", 40, fr);

            double[,] grille = sample2D(200, 200, scale, perlin(10000, TypeInterpolation.C2, 7, 3, 0.4, 10, 1.1, 0, 0.7));
            Console.WriteLine(minMax(grille));
            visu2D(grille,false).Save("D:/lab/Terrain/0perlinSimple.bmp");
            visu2D(erosionTerrainCroisee(ref r, ref grille,new double[] {0.2,1}, 80000, 0.5, 0.0002, 3,0.05, 0.05,2600,scale),true).Save("D:/lab/Terrain/Sediments.bmp");
            visu2D(grille,false).Save("D:/lab/Terrain/perlinErosion.bmp");
            
           grille = lissageMaritime(grille, 4,0.4);
            grille = lissage(grille, 1);
            visu2D(grille,false).Save("D:/lab/Terrain/perlinErosionLissee.bmp");
            exportCsv(grille, "D:/lab/Terrain/perlinErosionEau.csv");
            double[,] grillePente = sample2D(200, 200, scale, anglePente(grille, scale));
            visu2D(grillePente,false).Save("D:/lab/Terrain/perlinPente.bmp");
            double[,] ouv = ouvertureGrille(grille, scale, 0.0);
            visu2D(ouv,false).Save("D:/lab/Terrain/perlinOuv.bmp");
            Console.WriteLine(minMax(grille));
            Console.WriteLine(minMax(grillePente));
            Console.WriteLine(minMax(ouv));
            double[,] vege = couvertureVegetaleGrille(grille, ouv, grillePente, 0.05);
            visu2D(vege,false).Save("D:/lab/Terrain/vegetale.bmp");
            **/
        }

        public delegate double Application(double[] x);
        public delegate Bitmap Frame();
        public delegate void Del();

        //Fonctions pseudoAlea
        private static double pseudo_Alea(int N, int seed)
        {
            //Retourne un double pseudo aleatoire N->R
            N = N + seed*58900;
            N = (N << 13) ^ N;
            N = (N * (N * N * 15731 + 789221)) + 1376312589;
            return 1.0 - (N & 0x7fffffff) / 1073741824.0;
        }
        private static double pseudo_Alea_Rn(int[] V, int seed)
        {
            //Retourne un double pseudo aleatoire N^k->R
            int Taille = V.Length;
            double tmp = 0.0;
            for (int i = 0; i < Taille; i++)
            {
                tmp = (tmp * 850000.0);
                tmp = pseudo_Alea((int)tmp + V[i], seed);
            }
            return tmp;
        }
        private static double[] vect_pseudo_Alea_Rn(int[] V, int seed)
        {
            //Renvoie un vecteur de l'hypershere de dimention n, n etant la dimention de v
            int Taille = V.Length;
            double[] Retour = new double[Taille];
            double Norme;
            do
            {
                //Generer un vecteur dans l'hypercube
                for (int i = 0; i < Taille; i++)
                {
                    Retour[i] = pseudo_Alea_Rn(V, seed);
                    seed += 1;
                }
                Norme = normeVect(Retour);
            }
            while (Norme == 0 || Norme > 1.0);
            //Le vecteur est dans l'hyperboule, le normaliser
            multiplierVect(ref Retour, 1.0 / Norme);
            return Retour;
        }
        private static double[] vect_pseudo_Alea(int N,int Taille, int seed)
        {
            //Renvoie un vecteur de l'hypershere de dimention Taille
            double[] Retour = new double[Taille];
            double Norme;
            do
            {
                //Generer un vecteur dans l'hypercube
                for (int i = 0; i < Taille; i++)
                {
                    Retour[i] = pseudo_Alea(N, seed);
                    seed += 1;
                }
                Norme = normeVect(Retour);
            }
            while (Norme == 0 || Norme > 1.0);
            //Le vecteur est dans l'hyperboule, le normaliser
            multiplierVect(ref Retour, 1.0 / Norme);
            return Retour;
        }

        //Interpolation
        private static double fonctionsInterp(double x, TypeInterpolation Interp)
        {
            //Renvoie le coefficient d'interpolation k [0,1] en fonction du parcours x
            double result;
            switch (Interp)
            {
                case TypeInterpolation.Lineaire:
                    result = x;
                    break;
                case TypeInterpolation.Cosinus:
                    result = 0.5 * (1.0 - Math.Cos(x*Math.PI));
                    break;
                case TypeInterpolation.Hermite:
                    result = 3.0 * x * x - 2.0 * x * x * x;
                    break;
                case TypeInterpolation.C2:
                    result = 6.0 * x * x * x * x * x - 15.0 * x * x * x * x + 10.0 * x * x * x;
                    break;
                default:
                    result = 0.0;
                    break;
            }
            return result;
        }
        private static double interpolerEspace(double[] Position, int[] PosGrille, int N0, int seed, TypeInterpolation TI)
        {
            //Permet d'obtenir une interpolation pour une sous dim, a partir d'une selection des coordonees du point de ref pour les dim superieures
            if (N0 < 0)
            {
                //La valeur voulue pour un point de la grille est le produit scalaire entre le gradient a la position de la grille voulue, et la distance a cette position
                double[] dist = recentrer(PosGrille, Position);
                return produitScalaire(dist, vect_pseudo_Alea_Rn(PosGrille, seed));
            }
            else
            {
                //On choisit la dim N0 (2 choix) et on interpole les deux options entre elles
                double x = Position[N0];
                //projection sur la dim N0 et recuperation du coeff d'interpolation
                double k = fonctionsInterp(x - Math.Floor(x), TI);
                int[] P1 = copie(PosGrille);
                int[] P2 = copie(PosGrille);
                //fixer la valeur pour la dim N0
                P1[N0] = (int)Math.Floor(Position[N0]);
                P2[N0] = (int)Math.Floor(Position[N0]) + 1;
                //interpolation
                double b1 = interpolerEspace(Position, P1, N0 - 1, seed, TI);
                double b2 = interpolerEspace(Position, P2, N0 - 1, seed, TI);
                double inter = ((1.0 - k) * b1 + k * b2);
                return inter;
            }
        }
        private static double perlinSimple(ref double[] Position, int seed, TypeInterpolation TI)
        {
            //Obtient un bruit de perlin itere une fois et de frequence 1
            return Math.Min(Math.Max(interpolerEspace(Position, new int[Position.Length], Position.Length-1, seed, TI),-1.0),1.0);
        }
        public static double perlin(double[] Position, int seed, TypeInterpolation TI,int nbOctaves,double f0, double Attenuation,double Decalage,double puissance, int nbPlateaux, double k_plateaux)
        {
            //Renvoie un bruit de perlin a n dim pour l'interpolation specifiee, lissage par la puisssance "puissance"
            //nbPlateaux, nombre de plages de valeurs prises possibles. mettre à 0 pour toutes les valeurs possibles
            //k_plateaux est la pente des plateaux [0,1]
            double Resultat = 0.0;
            double Amplitude = 1.0;
            double f = f0;
            double[] shift = new double[Position.Length];
            double[] Pos = new double[Position.Length];
            double sommeAmp = 0;
            for (int i=0;i<nbOctaves;i++)
            {
                shift = vect_pseudo_Alea(i*452237+700849,Position.Length,seed);
                multiplierVect(ref shift, Decalage * pseudo_Alea(i*89746+6577,seed));
                Pos = copie(Position);
                shifter(ref Pos, ref shift);
                multiplierVect(ref Pos, f);
                Resultat += perlinSimple(ref Pos, seed, TI)*Amplitude;
                sommeAmp += Amplitude;
                Amplitude *= Attenuation;
                f *= 2;
            }
            Resultat = Resultat / sommeAmp;
            Resultat = Math.Sign(Resultat) * Math.Pow(Math.Abs(Resultat), puissance);
            if(nbPlateaux > 0)
            {
                double s = 0;
                if(k_plateaux != 0)
                {
                    s = k_plateaux * (1.0 / (float)nbPlateaux) * ((float)nbPlateaux * Resultat - Math.Floor((float)nbPlateaux * Resultat));
                }
                Resultat = (1.0 / (float)nbPlateaux) * Math.Floor((float)nbPlateaux * Resultat);
                Resultat += s;
            }
            return Resultat;
        }
        public static double perlin(double[] Position, int seed, TypeInterpolation TI, int nbOctaves, double f0, double Attenuation, double Decalage, double puissance,int nbPlateaux)
        {
            //Renvoie un bruit de perlin a n dim pour l'interpolation specifiee, lissage par la puisssance "puissance"
            //nbPlateaux, nombre de valeurs prises possibles
            return perlin(Position, seed, TI, nbOctaves, f0, Attenuation, Decalage, puissance, nbPlateaux, 0.0);
        }
        public static double perlin(double[] Position, int seed, TypeInterpolation TI, int nbOctaves, double f0, double Attenuation, double Decalage, double puissance)
        {
            //Renvoie un bruit de perlin a n dim pour l'interpolation specifiee, lissage par la puisssance "puissance"
            return perlin(Position, seed, TI, nbOctaves, f0, Attenuation, Decalage, puissance, 0, 0);
        }
        public static double perlin(double[] Position, int seed, TypeInterpolation TI, int nbOctaves, double f0, double Attenuation, double Decalage)
        {
            //Renvoie un bruit de perlin a n dim pour l'interpolation specifiee
            return perlin(Position, seed, TI, nbOctaves, f0, Attenuation, Decalage, 1.0, 0, 0);
        }
        public static double perlin(double[] Position, int seed, int nbOctaves, double f0, double Attenuation, double Decalage)
        {
            //Renvoie un bruit de perlin a n dim pour l'interpolation C2
            return perlin(Position, seed, TypeInterpolation.C2, nbOctaves, f0, Attenuation, Decalage, 1.0, 0, 0);
        }
        public static Application perlin(int seed, TypeInterpolation TI, int nbOctaves, double f0, double Attenuation, double Decalage, double puissance, int nbPlateaux, double k_plateaux)
        {
            return x => (1.0+perlin(x,seed,TI,nbOctaves,f0,Attenuation,Decalage,puissance,nbPlateaux,k_plateaux))/2.0;
        }

        //Fonctions vectorielles
        private static double produitScalaire(double[] V1, double[] V2)
        {
            //Produit scalaire de deux vecteurs
            int Taille = Math.Min(V1.Length, V2.Length);
            double Somme = 0.0;
            for (int i = 0; i < Taille; i++)
            {
                Somme += V1[i] * V2[i];
            }
            return Somme;
        }
        private static double normeVect(double[] V)
        {
            //Norme du vecteur v
            return Math.Sqrt(produitScalaire(V, V));
        }
        private static void multiplierVect(ref double[] V, double k)
        {
            //Multiplier le vecteur par k
            int Taille = V.Length;
            for (int i = 0; i < Taille; i++)
            {
                V[i] *= k;
            }
        }
        private static void shifter(ref double[] V1, ref double[] V2)
        {
            //renvoie le vecteur V1 = V1 +V2
            int Taille = Math.Min(V1.Length, V2.Length);
            for (int i = 0; i < Taille; i++)
            {
                V1[i] += V2[i];
            }
        }
        private static double[] recentrer(int[] V1, double[] V2)
        {
            //renvoie le vecteur V2-V1
            int Taille = Math.Min(V1.Length, V2.Length);
            double[] Resultat = new double[Taille];
            for (int i = 0; i < Taille; i++)
            {
                Resultat[i] = V2[i]-(double)V1[i];
            }
            return Resultat;
        }
        private static int[] copie(int[] V1)
        {
            //renvoie une copie de V1
            int Taille = V1.Length;
            int[] Resultat = new int[Taille];
            for (int i = 0; i < Taille; i++)
            {
                Resultat[i] = V1[i];
            }
            return Resultat;
        }
        private static double[] copie(double[] V1)
        {
            //renvoie une copie de V1
            int Taille = V1.Length;
            double[] Resultat = new double[Taille];
            for (int i = 0; i < Taille; i++)
            {
                Resultat[i] = V1[i];
            }
            return Resultat;
        }

        //Fonctions affichage
        private static double[,] sample2D(int x, int y, double scale, Application f)
        {
            double[,] retour = new double[x, y];
            double[] x_;
            for (int i = 0; i < x; i++)
            {
                for (int j = 0; j < y; j++)
                {
                    x_ = new double[] { i * scale, j * scale };
                    double y_ = f(x_);

                    retour[i, j] = y_;
                }
            }
            return retour;
        }
        private static Tuple<double,double> minMax(double[,]grille)
        {
            int x = grille.GetLength(0);
            int y = grille.GetLength(1);
            double max = double.NegativeInfinity;
            double min = double.PositiveInfinity;
            for (int i = 0; i < x; i++)
            {
                for (int j = 0; j < y; j++)
                {
                    if (grille[i, j] < min)
                    {
                        min = grille[i, j];
                    }
                    if (grille[i, j] > max)
                    {
                        max = grille[i, j];
                    }
                }
            }
            return new Tuple<double, double>(max, min);
        }
        private static void Normer(double[,] grille)
        {
            int x = grille.GetLength(0);
            int y = grille.GetLength(1);
            Tuple<double, double> minmax = minMax(grille);
            double max = minmax.Item1;
            double min = minmax.Item2;
            for (int i = 0; i < x; i++)
            {
                for (int j = 0; j < y; j++)
                {
                    grille[i, j] = (grille[i, j] - min) / (max - min);
                }
            }

        }
        private static Bitmap visu2D( double[,] grille,bool discret)
        {
            int x = grille.GetLength(0);
            int y = grille.GetLength(1);
            Bitmap retour = new Bitmap(x, y);
            double[,] gN = (double[,])grille.Clone();
            Normer(gN);
            for(int i=0;i<x;i++)
            {
                for(int j=0;j<y;j++)
                {

                    double y_ = gN[i, j];
                    Color c;
                    double hue = y_ ;
                    int val;
                    val = Math.Min(255, Math.Max(0, (int)(hue * 255.0)));
                    if(discret)
                    {
                        if(val == 0)
                        {
                            c = Color.FromArgb(0, 0, 0);
                        }
                        else
                        {
                            c = Color.FromArgb(255, 255, 255);
                        }
                        
                    }
                    else
                    {
                        c = Color.FromArgb(val, val, val);
                    }
                    retour.SetPixel(i, j, c);
                }
            }
            return retour;
        }
        private static Bitmap visu2D(int[,] grille)
        {
            int x = grille.GetLength(0);
            int y = grille.GetLength(1);
            Bitmap retour = new Bitmap(x, y);
            for (int i = 0; i < x; i++)
            {
                for (int j = 0; j < y; j++)
                {

                    Color c;
                    if(grille[i,j] ==0)
                    {
                        c = Color.Gray;
                    }
                    else if (grille[i, j] == -1)
                    {
                        c = Color.Blue;
                    }
                    else if (grille[i, j] == 1)
                    {
                        c = Color.Red;
                    }
                    else
                    {
                        c = Color.Black;
                    }
                    retour.SetPixel(i, j, c);
                }
            }
            return retour;
        }
        private static void exportCsv(double[,] grille,string ExportPath)
        {
            int x = grille.GetLength(0);
            int y = grille.GetLength(1);
            using (System.IO.StreamWriter file =
           new System.IO.StreamWriter(ExportPath))
            {
                for (int i = 0; i < x; i++)
                {
                    String ligne = "";
                    for (int j = 0; j < y; j++)
                    {

                        double y_ = grille[i, j];
                        ligne += y_ + ";";
                    }
                    file.WriteLine(ligne);
                }
            }
        }
        private static void creerGIF(string Path,int nbFrame , Frame frame)
        {

            System.Windows.Media.Imaging.GifBitmapEncoder gEnc = new GifBitmapEncoder();
            for (int i = 0; i < nbFrame; i++)
            {
                var bmp = frame().GetHbitmap();
                var src = System.Windows.Interop.Imaging.CreateBitmapSourceFromHBitmap(
                        bmp,
                        IntPtr.Zero,
                        System.Windows.Int32Rect.Empty,
                        BitmapSizeOptions.FromEmptyOptions());
                        gEnc.Frames.Add(BitmapFrame.Create(src));

            }
            using (System.IO.FileStream fs = new System.IO.FileStream(Path, System.IO.FileMode.Create))
            {
                gEnc.Save(fs);
            }
        }


        //Erosion
        private static double[] gradient(double[,] grille, double[] x, double scale)
        {
            int X = grille.GetLength(0);
            int Y = grille.GetLength(1);
            int i = Math.Min(X - 1, Math.Max(0, (int)Math.Floor(x[0] / scale)));
            int j = Math.Min(Y - 1, Math.Max(0, (int)Math.Floor(x[1] / scale)));
            double xp = grille[Math.Min(X - 1, i + 1), j];
            double xm = grille[Math.Max(0, i - 1), j];
            double yp = grille[i, Math.Min(Y - 1, j + 1)];
            double ym = grille[i, Math.Max(0, j - 1)];
            double dx1 = ( xp-xm ) / (2.0 * scale);
            double dy1 = ( yp- ym) / (2.0 * scale);
            return new double[] { -(dx1) , -(dy1) };
        }
        private static double[] gradient(double[,] grille, double[] x, double scale,ref double gradLen)
        {
            double[] gr = gradient(grille, x, scale);
            double dx = gr[0];
            double dy = gr[1];
            double len = Math.Sqrt(dx*dx+dy*dy);
            gradLen = len;
            if(len ==0)
            {
                return new double[] { 0, 0 };
            }
            return new double[] { dx/len, dy/len };
        }
        private static void afficherVect(double[] x)
        {
            Console.WriteLine("[" + x[0] + " " + x[1] + "]");
        }
        private static double[] incrementer(double[] x, double[] dx, double delta)
        {
            return new double[] { x[0] + dx[0] * delta, x[1] + dx[1] * delta };
        }
        private static double[] lerp(double[] x0, double[] x1, double k)
        {
            double xn = x0[0] * k + x1[0] * (1.0 - k);
            double yn = x0[1] * k + x1[1] * (1.0 - k);
            double len = Math.Sqrt(xn * xn + yn * yn);
            if(len ==0)
            {
                return new double[] {0 , 0}; ;
            }
            return new double[] {xn/len , yn/len };
        }
        private static double norme(double[]x)
        {
            return Math.Sqrt(Math.Pow(x[0], 2) + Math.Pow(x[1], 2));
        }
        private static double[] biaiser(ref Random r, double[] x, double biais)
        {
            double len = 0;
            double[] retour = new double[2];
            while (len == 0)
            {
                double teta = r.NextDouble() * Math.PI * 2.0;
                double rayon = biais * r.NextDouble();
                double bx = rayon * Math.Cos(teta);
                double by = rayon * Math.Sin(teta);
                retour = new double[] { x[0] + bx, x[1] +  by};
                len = norme(retour);
                if (len != 0)
                {
                    retour[0] /= len;
                    retour[1] /= len;
                }
            }
            return retour;
        }
        private static double[,] erosionTerrainCroisee(ref Random r, ref double[,] grille, double[] vent, int nbErosions, double capaciteMax, double evaporation, double erosion, double biais, double inertie, double ratioEauVent, double scale)
        {
            int X = grille.GetLength(0);
            int Y = grille.GetLength(1);
            double[,] sediments = new double[X, Y];
            for (int i = 0; i < X; i++)
            {
                for (int j = 0; j < Y; j++)
                {
                    sediments[i, j] = 0;
                }
            }
            for (int i = 0; i < nbErosions; i++)
            {
                erosionVent(ref r, ref grille,ref sediments, vent, capaciteMax / (ratioEauVent), evaporation/(ratioEauVent), erosion/ratioEauVent, biais, scale);
                erosionGoutte(ref r, ref grille,ref sediments, capaciteMax, evaporation, erosion, scale, inertie,biais);
                if(i%1000==0)
                {
                    Console.WriteLine(i);
                }
            }
            return sediments;
        }
        private static void erosionGoutte(ref Random r, ref double[,] grille,ref double[,] sediments, double capaciteMax, double evaporation, double erosion, double scale, double inertie, double biais)
        {
            int X = grille.GetLength(0);
            int Y = grille.GetLength(1);
            double[] x = new double[] { X * r.NextDouble() * scale, Y * r.NextDouble() * scale };
            double[] pas = new double[] { 0, 0 };
            double qteSediment = 0;
            double capacite = r.NextDouble() * capaciteMax;
            double deltaH = 0;
            double saut = 0.25;
            while (capacite>0 || qteSediment >0)
            {
                
                int ir = (int)Math.Floor(x[0] / scale);
                int jr = (int)Math.Floor(x[1] / scale);
                if(ir<0 || jr<0 || ir>(X-1)||jr>(Y-1))
                {
                    break;
                }
                int i = Math.Min(X - 1, Math.Max(0, ir));
                int j = Math.Min(Y - 1, Math.Max(0, jr));
                double gradLen=0;
                double[] gradient_ = gradient(grille, x, scale,ref gradLen);
                pas = lerp(pas,gradient_ , inertie);
                double[] xnext;
                int inext;
                int jnext;
                
                pas = biaiser(ref r, pas, biais);
                xnext = incrementer(x, pas, saut * scale);
                inext = Math.Min(X - 1, Math.Max(0, (int)Math.Floor(xnext[0] / scale)));
                jnext = Math.Min(Y - 1, Math.Max(0, (int)Math.Floor(xnext[1] / scale)));
                deltaH = grille[i, j] - grille[inext, jnext];
                if (capacite>qteSediment)
                {
                    
                    double qteErosion = Math.Min(Math.Max(0, deltaH), erosion*gradLen);
                    double transfert = Math.Min(qteErosion, capacite - qteSediment);
                    qteSediment += transfert;
                    sediments[i, j] -= Math.Min(transfert, sediments[i, j]);
                    grille[i,j] -= transfert;
                }
                else
                {
                    double transfert = Math.Min((qteSediment - capacite), 0.001+Math.Max(0,-deltaH));
                    qteSediment -=transfert;
                    sediments[i, j] += transfert;
                    grille[i,j] += transfert;
                }
                capacite = Math.Max(0, capacite - evaporation);
                x = xnext;
                //Console.WriteLine(capacite + "/" + qteSediment);

            }
        }
        private static void erosionVent(ref Random r, ref double[,] grille, ref double[,] sediments,double[] vent, double capaciteMax, double evaporation, double erosion,double biais, double scale)
        {
            int X = grille.GetLength(0);
            int Y = grille.GetLength(1);
            double[] x = new double[] { X * r.NextDouble() * scale, Y * r.NextDouble() * scale };
            multiplierVect(ref vent, 1.0 / norme(vent));
            vent = biaiser(ref r, vent, biais);
            double qteSediment = 0;
            double capacite = r.NextDouble() * capaciteMax;
            while (capacite > 0 || qteSediment>0)
            {
                int ir = (int)Math.Floor(x[0] / scale);
                int jr = (int)Math.Floor(x[1] / scale);
                if (ir < 0 || jr < 0 || ir > (X - 1) || jr > (Y - 1))
                {
                    break;
                }
                int i = Math.Min(X - 1, Math.Max(0, ir));
                int j = Math.Min(Y - 1, Math.Max(0, jr));
                double[] xnext;
                int inext;
                int jnext;
                double saut = 0.25;
                xnext = incrementer(x, vent, saut * scale);
                inext = Math.Min(X - 1, Math.Max(0, (int)Math.Floor(xnext[0] / scale)));
                jnext = Math.Min(Y - 1, Math.Max(0, (int)Math.Floor(xnext[1] / scale)));
                double deltaH = grille[i, j] - grille[inext, jnext];
                if (capacite > qteSediment)
                {
                    double qteErosion = Math.Min(Math.Max(0,deltaH), erosion);
                    double transfert = Math.Min(qteErosion, capacite - qteSediment);
                    qteSediment += transfert;
                    sediments[i, j] -= Math.Min(transfert, sediments[i, j]);
                    grille[i, j] -= transfert;
                }
                else
                {
                    double transfert = (qteSediment - capacite);
                    qteSediment = capacite;
                    sediments[i, j] += transfert;
                    grille[i, j] += transfert;
                }
                capacite = Math.Max(0, capacite - evaporation);
                x = xnext;
            }
        }

        //Pente
        public static double[] normale(double[] grad)
        {
            //renvoie le vecteur normal (si f va de r^n dans r, la normale est un vecteur de r^n+1) conaissant le gradient
            //le vecteur normal est normal a toutes les derivees partielles du vecteur M(x,y,...,f(x,y,..))
            int Taille = grad.Length;
            double[] retour = new double[Taille + 1];
            retour[Taille] = 1.0;
            for (int i = 0; i < Taille; i++)
            {
                retour[i] = -grad[i];
            }
            double norme = normeVect(retour);
            multiplierVect(ref retour, 1.0 / norme);
            return retour;
        }
        public static double anglePente(double[] normale)
        {
            //Renvoie l'inclinaison en rad de l'hyperplan tangent à la fonction conaissant la normale
            double cos = normale[normale.Length - 1];//produit scalaire avec la derniere dim
            double alpha = Math.Acos(cos);
            return alpha;
        }
        public static double anglePente(double[] x, double[,] grille, double scale)
        {
            //Renvoie l'inclinaison en rad de l'hyperplan tangent à la fonction conaissant la fonction
            return anglePente(normale(gradient(grille,x,scale)));
        }
        public static Application anglePente( double[,] grille, double scale)
        {
            return (x) => Math.Abs(anglePente(x,grille, scale))/Math.PI;
        }

        //Exposition
        public static double occlusionLocale(double[,] grille, int i0, int j0, int i,int j, double scale)
        {
            double dx = Math.Sqrt(Math.Pow(i - i0, 2.0) + Math.Pow(j - j0, 2.0)) * scale;
            double dy = Math.Max(0, grille[i, j] - grille[i0, j0]);
            if (dx != 0)
            {
                return Math.Atan(dy / dx);
            }
            else
            {
                return 0;
            }
        }
        public static double ouverture(double[,] grille, int i0, int j0, double scale, double angleSoleil)
        {
            return 1.0 - (occlusionMax(grille, i0, j0, scale, angleSoleil, 1) + occlusionMax(grille, i0, j0, scale, angleSoleil, -1))/Math.PI; 
        }
        public static double occlusionMax(double[,] grille, int i0, int j0, double scale,double angleSoleil,double signe)
        {
            int X = grille.GetLength(0);
            int Y = grille.GetLength(1);
            double k = 0;
            int ir = i0;
            int jr = j0;
            double occlusionMx = 0;
            while(true)
            {
                ir = (int)(i0 + signe*k * Math.Cos(angleSoleil));
                jr = (int)(j0 + signe*k * Math.Sin(angleSoleil));
                if (ir < 0 || jr < 0 || ir > (X - 1) || jr > (Y - 1))
                {
                    break;
                }
                k += 0.5;
                double occlusion = occlusionLocale(grille, i0, j0, ir, jr, scale);
                if(occlusion>occlusionMx)
                {
                    occlusionMx = occlusion;
                }
            }
            return occlusionMx;
        }
        public static double[,] ouvertureGrille(double[,] grille, double scale, double angleSoleil)
        {
            int X = grille.GetLength(0);
            int Y = grille.GetLength(1);
            double[,] retour = new double[X, Y];
            for (int i = 0; i < X; i++)
            {
                for (int j = 0; j < Y; j++)
                {
                    double ouv = ouverture(grille, i, j, scale, angleSoleil);
                    retour[i, j] = ouv;
                }
            }
            return retour;
        }

        //vegetation
        public static double expInverse(double x, double lambda)
        {
            return  Math.Exp(-x * lambda);
        }
        public static double couvertureVegetale(double altitude, double ouverture, double pente, double lambdaAltitude, double lambdaCouverture, double lambdaPente)
        {
            double compoAlt = expInverse(altitude, lambdaAltitude);
            double compoOuv = expInverse(1.0 - ouverture, lambdaCouverture);
            double compoPente = expInverse(pente, lambdaPente); ;
            return compoAlt * compoOuv * compoPente;
        }
        public static double[,] couvertureVegetaleGrille(double[,] altitude, double[,] ouverture, double[,] pente, double lambdaAltitude, double lambdaCouverture, double lambdaPente)
        {
            int X = altitude.GetLength(0);
            int Y = altitude.GetLength(1);
            double[,] retour = new double[X, Y];
            for (int i = 0; i < X; i++)
            {
                for (int j = 0; j < Y; j++)
                {
                    double couv = couvertureVegetale(altitude[i,j],ouverture[i,j],pente[i,j],lambdaAltitude,lambdaCouverture,lambdaPente);
                    retour[i, j] = couv;
                }
            }
            return retour;
        }
        public static double[,] couvertureVegetaleGrille(double[,] altitude, double[,] ouverture, double[,] pente,double attenationMax)
        {
            double maxAlt = minMax(altitude).Item2;
            double couvMax = 1.0-minMax(ouverture).Item1;
            double penteMax = minMax(pente).Item2;
            double sigmaAlt = -Math.Log(attenationMax) / (3*maxAlt);
            double sigmaCouv = -Math.Log(attenationMax) / (3 * couvMax);
            double sigmaPente = -Math.Log(attenationMax) / (3 * penteMax);
            return couvertureVegetaleGrille(altitude, ouverture, pente, sigmaAlt, sigmaCouv, sigmaPente);
        }
        //Lissage
        public static double[,] lissageMaritime(double[,]grille , int rLissage,double zMax)
        {
            int X = grille.GetLength(0);
            int Y = grille.GetLength(1);
            double[,] retour = new double[X, Y];
            for (int i = 0; i < X; i++)
            {
                for (int j = 0; j < Y; j++)
                {
                    double somme = 0;
                    double n = 0;
                    double local = grille[i, j];
                    if (local <= zMax)
                    {
                        for (int k = -rLissage; k <= rLissage; k++)
                        {
                            for (int l = -rLissage; l <= rLissage; l++)
                            {
                                double r = Math.Sqrt((double)(l * l) + (double)(k * k));
                                double coeff;
                                coeff = 1.0 / (1.0 + r);
                                int inext = Math.Min(X - 1, Math.Max(0, i + k));
                                int jnext = Math.Min(Y - 1, Math.Max(0, j + l));
                                double val = grille[inext, jnext];
                                somme += coeff * val;
                                n += coeff;
                            }
                        }
                        retour[i, j] = somme / n;
                    }
                    else
                    {
                        retour[i, j] = local;
                    }
                }
            }
            return retour;
        }
        public static double[,] lissage(double[,] grille, int rLissage)
        {
            Tuple<double, double> minmax = minMax(grille);
            double max = minmax.Item1;
            return lissageMaritime(grille, rLissage, max+1);
        }

        //Tectonique
        public static Application Cos()
        {
            return (x) => Math.Cos(x[0]*Math.PI*2);
        }
        public static Application Cassure(int seed)
        {
            return (x) => {
                double k = x[0];
                double partDec = k - Math.Floor(k);
                if(partDec<=0.5)
                {
                    return partDec*2;
                }
                else
                {
                    double part = (partDec - 0.5) / 0.5;
                    
                    double A0 = (1.0+pseudo_Alea((int)Math.Floor(k), seed))/2.0;
                    return (1.0 - part) * A0;
                }
            };
        }
        public static Application Triangle()
        {
            return (x) => {
                double k = x[0];
                double partDec = k - Math.Floor(k);
                if (partDec <= 0.5)
                {
                    return partDec * 2;
                }
                else
                {
                    double part = (partDec - 0.5) / 0.5;
                    return (1.0 - part);
                }
            };
        }
        public static Application Ondelettes(double f, double[] centre)
        {
            return (x) =>
            {
                shifter(ref x, ref centre);
                return Math.Cos(normeVect(x) * f * Math.PI * 2);
                };
        }
        public static double[,] froisser(double[,] grille,double[] direction,double[] decalage, double SigmaAttenuationLaterale, double FrequenceColineaire, Application MotifColineaire, Application BruitMultiplicatif,Application MotifNormal,double scale,double amplitude,double AmplitudeMotifNormal,double frequence_decalageFront)
        {
            int X = grille.GetLength(0);
            int Y = grille.GetLength(1);
            double[,] retour = new double[X, Y];
            multiplierVect(ref direction, 1.0 / norme(direction));
            double[] x_;
            for (int i = 0; i < X; i++)
            {
                for (int j = 0; j < Y; j++)
                {
                    x_ = new double[] { i * scale-decalage[0], j * scale-decalage[1] };
                    double proj = x_[0] * direction[0] + x_[1] * direction[1];
                    double[] ecartement = new double[] { x_[0] - proj * direction[0], x_[0] - proj * direction[0] };
                    double dist = norme(ecartement);
                    double posFront = proj - MotifNormal(new double[] {dist*Math.Sign(ecartement[0])-proj*frequence_decalageFront})*AmplitudeMotifNormal;
                    double y_ = MotifColineaire(new double[] { posFront*FrequenceColineaire }) * Math.Exp(-Math.Pow(dist / SigmaAttenuationLaterale, 2)) * BruitMultiplicatif(x_);
                    retour[i, j] =grille[i,j]+  y_*amplitude;
                }
            }
            Normer(retour);
            return retour;
        }

        //Spatial
        public static double HauteurCratere(double rayon, double R0, double RaideurCassure,double HauteurCassure, double RaideurDepression, double AmplitudeDepression, double RaideurRebond, double H0)
        {
            double depress = -AmplitudeDepression / (1.0 + Math.Exp(RaideurDepression * (rayon - R0)));
            double cass = Math.Exp(-Math.Abs(rayon - R0)*RaideurCassure)*(HauteurCassure+AmplitudeDepression/2.0);
            double h0 = -AmplitudeDepression / (1.0 + Math.Exp(RaideurDepression * (-R0))) + Math.Exp(-Math.Abs(-R0)) * (HauteurCassure + AmplitudeDepression / 2.0);
            double reb = Math.Exp(-RaideurRebond * Math.Pow(rayon, 2)) * (H0 - h0);
            return depress + cass + reb;
        }
        public static double[,] Impact(double[,] grille, double scale, double[] position, Application BruitRayon,double AmplitudeBruit, double R0, double RaideurCassure, double HauteurCassure, double RaideurDepression, double AmplitudeDepression, double RaideurRebond, double H0,double AttenuationLissage, double RayonLissageRelatif)
        {
            int X = grille.GetLength(0);
            int Y = grille.GetLength(1);
            double[,] retour = new double[X, Y];
            int RayonLissage = (int)(RayonLissageRelatif * R0 / scale);
            double[,] lisse = lissage(grille,RayonLissage );
            
            for (int i = 0; i < X; i++)
            {
                for (int j = 0; j < Y; j++)
                {
                    double[] x_ = new double[] { i * scale, j * scale};
                    double r = Math.Sqrt(Math.Pow(x_[0] - position[0], 2) + Math.Pow(x_[1] - position[1], 2)) + BruitRayon(x_)*AmplitudeBruit;
                    double h = HauteurCratere(r, R0, RaideurCassure, HauteurCassure, RaideurDepression, AmplitudeDepression, RaideurRebond, H0);
                    double k_lissage = Math.Exp(-AttenuationLissage * r/R0);
                    retour[i, j] = (1.0-k_lissage)*grille[i, j] + k_lissage*lisse[i,j] + h;
                }
            }
            return retour;
        }
        public static double[,] ImpactAleatoire(double[,] grille, double scale,ref Random r, Application BruitRayon,double AmplitudeBruitMax, double R0Max, double RaideurCassureMax, double HauteurCassureMax, double RaideurDepressionMax, double AmplitudeDepressionMax, double RaideurRebondMax, double H0Max, double AttenuationLissageMax, double RayonLissageRelatifMax)
        {
            int X = grille.GetLength(0);
            int Y = grille.GetLength(1);
            double[] position = new double[] {X*scale*r.NextDouble(), Y * scale * r.NextDouble() };
            return Impact(grille, scale, position, BruitRayon,AmplitudeBruitMax*r.NextDouble(), R0Max * r.NextDouble(), RaideurCassureMax * r.NextDouble(), HauteurCassureMax * r.NextDouble(),  RaideurDepressionMax * r.NextDouble(), AmplitudeDepressionMax * r.NextDouble(), RaideurRebondMax * r.NextDouble(), H0Max * r.NextDouble(), AttenuationLissageMax*r.NextDouble(), RayonLissageRelatifMax*r.NextDouble());
        }
        public static double[,] Bombarder(double[,] grille, int nbImpacts, double scale, ref Random r, Application BruitRayon, double AmplitudeBruitMax, double R0Max, double RaideurCassureMax, double HauteurCassureMax, double RaideurDepressionMax, double AmplitudeDepressionMax, double RaideurRebondMax, double H0Max, double AttenuationLissageMax, double RayonLissageRelatifMax)
        {
            double[,] resultat = (double[,])grille.Clone();
            for(int i=0;i<nbImpacts;i++)
            {
                Console.WriteLine(i);
                resultat = ImpactAleatoire(resultat, scale, ref r, BruitRayon, AmplitudeBruitMax, R0Max, RaideurCassureMax, HauteurCassureMax, RaideurDepressionMax, AmplitudeDepressionMax, RaideurRebondMax, H0Max,AttenuationLissageMax,RayonLissageRelatifMax);
                visu2D(resultat, false).Save("D:/lab/Terrain/bombe" + i + ".bmp");
            }
            return resultat;
        }
        //Lacs
#pragma warning disable CS0659 // Le type se substitue à Object.Equals(object o) mais pas à Object.GetHashCode()
        public class Case
        {
            public int x;
            public int y;

            public Case(int x, int y)
            {
                this.x = x;
                this.y = y;
            }
            public Case(Case c)
            {
                this.x = c.x;
                this.y = c.y;
            }
            public override bool Equals(object obj)
            {
                var @case = obj as Case;
                return @case != null &&
                       x == @case.x &&
                       y == @case.y;
            }
            public bool bordure(int X, int Y)
            {
                return x == 0 || y == 0 || x == X - 1 || y == Y - 1;
            }
            
        }
        public static bool construireLac(ref List<Case> cases, ref double[,] altitudes, double zmax, int x, int y,int pr )
        {
            //if(pr %100 ==0)
            //Console.WriteLine(pr);
            int X = altitudes.GetLength(0);
            int Y = altitudes.GetLength(1);
            Case c = new Case(x, y);
            if(!cases.Contains(c))
            {
                if (altitudes[c.x,c.y]<zmax)
                {
                    cases.Add(c);
                    bool fuite = c.bordure(X, Y);
                    if (!fuite)
                    {
                        fuite = fuite || construireLac(ref cases, ref altitudes, zmax, x + 1, y,pr+1);
                    }
                    if (!fuite)
                    {
                        fuite = fuite || construireLac(ref cases, ref altitudes, zmax, x - 1, y,pr+1);
                    }
                    if(!fuite)
                    {
                        fuite = fuite || construireLac(ref cases, ref altitudes, zmax, x, y + 1,pr+1);
                    }
                    if (!fuite)
                    {
                        fuite = fuite || construireLac(ref cases, ref altitudes, zmax, x, y - 1,pr+1);
                    }
                    return fuite;
                }
                else
                {
                    return false;
                }

            }
            else
            {
                return false;
            }
        }
        public static void trouverAltitudeMaxLac(ref double[,] altitudes, ref double[,] altitudesLacs,double profMin, int x, int y,double deltaZ)
        {
            List<Case> cases = new List<Case>();
            List<Case> casesTest = new List<Case>();
            double z0 = altitudes[x, y] + profMin;
            double z = z0;
            while(!construireLac(ref casesTest,ref altitudes,z,x,y,0))
            {
                cases = casesTest;
                casesTest = new List<Case>();
                z += deltaZ;
            }
            if (z != z0)
            {
                foreach (Case c in cases)
                {
                    altitudesLacs[x, y] = z;
                }
            }
        }
        public static double[,] altitudesLacs_(double[,] altitudes, double profMin, double deltaZ)
        {
            int X = altitudes.GetLength(0);
            int Y = altitudes.GetLength(1);
            double[,] lacs = new double[X, Y];
            for (int i = 0; i < X; i++)
            {
                for (int j = 0; j < Y; j++)
                {
                    lacs[i, j] = 0;
                }
            }
            for (int i = 0; i < X; i++)
            {
                Console.WriteLine(i);
                for (int j = 0; j < Y; j++)
                {
                    if (lacs[i, j] == 0)
                    {
                        trouverAltitudeMaxLac(ref altitudes, ref lacs, profMin, i, j, deltaZ);
                    }
                }
            }
            return lacs;

        }
        public static double[,] altitudesLacs(double[,] altitudes, double profMin, double deltaZ)
        {
            int N = 100 * altitudes.Length;
            double[,] lacs = new double[0, 0];
            Del BigRecursion = () => {
                lacs = altitudesLacs_(altitudes, profMin, deltaZ);
            };
            Thread thread = new Thread(new ThreadStart(BigRecursion), N);
            thread.Start();
            thread.Join();
            return lacs;
        }
    }
}
