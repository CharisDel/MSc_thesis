#include <iostream> // For input and output operations
#include <vector>   // For using std::vector
#include <cmath>    // For mathematical functions like log
#include <TMath.h>  // For ROOT's mathematical functions (if used)
#include <TRandom3.h>
#include <TGraph.h>
#include "TLatex.h"
#include <TCanvas.h>
#include <TLine.h>

double calc_neglnL(double Nobs, double ns, double mu_s, double ntop, double mu_Top, double nWZ, double mu_WZ, double nZ2, double mu_Z2, double nZ1, double mu_Z1, double nZ0, double mu_Z0, double nWW, double mu_WW, double N_bkg)
{
    double lnLi;
    lnLi = +Nobs * log(mu_s * ns + mu_Top * ntop + mu_WZ * nWZ + mu_Z2 * nZ2 + mu_Z1 * nZ1 + mu_Z0 * nZ0 + mu_WW * nWW + N_bkg) - mu_s * ns - mu_Top * ntop - mu_WZ * nWZ - mu_Z2 * nZ2 - mu_Z1 * nZ1 - mu_Z0 * nZ0 - mu_WW * nWW - N_bkg;
    return -lnLi;
}

double calc_min(const std::vector<double> &ar, long int n)
{
    bool min = 0;
    double dummy = 0;
    dummy = ar[0];
    int dummy2;
    for (int i = 0; i <= n; i++)
    {
        if (dummy >= ar[i])
        {
            dummy = ar[i];
            dummy2 = i;
        }
    }
    return dummy2;
}

void Simultane_new()
{
    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    // different for every area
    long double lnL_signal1;
    long double lnL_signal2;
    long double lnL_top1;
    long double lnL_top2;
    long double lnL_3lCR1;
    long double lnL_3lCR2;
    long double lnL_Zjets21;
    long double lnL_Zjets22;
    long double lnL_Zjets11;
    long double lnL_Zjets12;
    long double lnL_Zjets01;
    long double lnL_Zjets02;
    long double lnL_WW1;
    long double lnL_WW2;
    long double lnL_all;

    double mu_s;
    double mu_top;
    double mu_WZ;
    double mu_Zjets2;
    double mu_Zjets1;
    double mu_Zjets0;
    double mu_WW;
    int min = 0;
    long int i = 1;
    int i_best = 0;

    double min_lnL = +9999;
    double Min_lnL = 9999;
    double Min_LnL = 9999;
    std::vector<double> lnL_tot;
    std::vector<double> lnL_tot1(100);
    std::vector<double> mu_s_dummy;
    std::vector<double> mu_s_dummy1(100);
    std::vector<double> mu_top_dummy;
    std::vector<double> mu_top_dummy1(100);
    std::vector<double> mu_WZ_dummy;
    std::vector<double> mu_WZ_dummy1(100);
    std::vector<double> mu_Zjets2_dummy;
    std::vector<double> mu_Zjets2_dummy1(100);
    std::vector<double> mu_Zjets1_dummy;
    std::vector<double> mu_Zjets1_dummy1(100);
    std::vector<double> mu_Zjets0_dummy;
    std::vector<double> mu_Zjets0_dummy1(100);
    std::vector<double> mu_WW_dummy;
    std::vector<double> mu_WW_dummy1(100);

    for (mu_top = 1.00; mu_top <= 1.02; mu_top += 0.002) // boundaries for every scaling factor are set +-5σ from its nominal value that has been found
    {
        for (mu_WZ = 1.00; mu_WZ <= 1.02; mu_WZ += 0.002)
        {
            for (mu_Zjets2 = 1.06; mu_Zjets2 <= 1.08; mu_Zjets2 += 0.002)
            {
                for (mu_Zjets1 = 1.31; mu_Zjets1 <= 1.33; mu_Zjets1 += 0.002)
                {
                    for (mu_Zjets0 = 1.34; mu_Zjets0 <= 1.36; mu_Zjets0 += 0.002)
                    {
                        for (mu_WW = 1.4; mu_WW <= 1.50; mu_WW += 0.002)
                        {
                            for (mu_s = 0.99; mu_s <= 1.01; mu_s += 0.001)
                            {
                                // obs  signal         top             WZ             Z2               Z1                Z0                 WW            Other
                                lnL_signal1 = calc_neglnL(1295, 667.74, mu_s, 78.38, mu_top, 339.55, mu_WZ, 1.22, mu_Zjets2, 39.51, mu_Zjets1, 72.03, mu_Zjets0, 18.52, mu_WW, 21.65);

                                lnL_signal2 = calc_neglnL(1426, 939.55, mu_s, 22.94, mu_top, 348.29, mu_WZ, 0.22, mu_Zjets2, 3.14, mu_Zjets1, 35.89, mu_Zjets0, 19.99, mu_WW, 24.32);

                                lnL_top1 = calc_neglnL(2867, 0.0, mu_s, 2860.72, mu_top, 0.38, mu_WZ, 0.06, mu_Zjets2, 0.00, mu_Zjets1, 0.0, mu_Zjets0, 6.07, mu_WW, 0.92);

                                lnL_top2 = calc_neglnL(2869, 0.0, mu_s, 2792.99, mu_top, 0.48, mu_WZ, 0.10, mu_Zjets2, 0.12, mu_Zjets1, 0.0, mu_Zjets0, 6.96, mu_WW, 3.64);

                                lnL_3lCR1 = calc_neglnL(1229, 0.86, mu_s, 34.00, mu_top, 1080.97, mu_WZ, 9.18, mu_Zjets2, 20.52, mu_Zjets1, 2.45, mu_Zjets0, 0.37, mu_WW, 51.30);

                                lnL_3lCR2 = calc_neglnL(1180, 1.42, mu_s, 22.97, mu_top, 1018.03, mu_WZ, 6.53, mu_Zjets2, 13.96, mu_Zjets1, 33.29, mu_Zjets0, 0.27, mu_WW, 62.80);

                                lnL_Zjets21 = calc_neglnL(2871, 83.64, mu_s, 152.35, mu_top, 158.64, mu_WZ, 2296.70, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 12.44, mu_WW, 24.86);

                                lnL_Zjets22 = calc_neglnL(2898, 78.29, mu_s, 91.75, mu_top, 150.05, mu_WZ, 2362.55, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 5.95, mu_WW, 33.03);

                                lnL_Zjets11 = calc_neglnL(2836, 86.51, mu_s, 144.00, mu_top, 144.37, mu_WZ, 0.0, mu_Zjets2, 1790.57, mu_Zjets1, 0.0, mu_Zjets0, 34.78, mu_WW, 18.44);

                                lnL_Zjets12 = calc_neglnL(2897, 58.53, mu_s, 52.53, mu_top, 103.86, mu_WZ, 0.0, mu_Zjets2, 2007.90, mu_Zjets1, 0.0, mu_Zjets0, 11.20, mu_WW, 27.66);

                                lnL_Zjets01 = calc_neglnL(1355, 36.96, mu_s, 63.07, mu_top, 40.47, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 855.74, mu_Zjets0, 19.71, mu_WW, 7.57);

                                lnL_Zjets02 = calc_neglnL(1336, 60.77, mu_s, 27.94, mu_top, 38.90, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 874.37, mu_Zjets0, 27.96, mu_WW, 4.84);

                                lnL_WW1 = calc_neglnL(1445, -0.08, mu_s, 1067.54, mu_top, 8.99, mu_WZ, 0.19, mu_Zjets2, 2.29, mu_Zjets1, 0.29, mu_Zjets0, 206.19, mu_WW, 28.21);

                                lnL_WW2 = calc_neglnL(1458, 0.49, mu_s, 767.31, mu_top, 13.30, mu_WZ, 0.21, mu_Zjets2, -0.07, mu_Zjets1, 0.35, mu_Zjets0, 425.87, mu_WW, 68.67);

                                lnL_all = lnL_signal1 + lnL_signal2 + lnL_top1 + lnL_top2 + lnL_3lCR1 + lnL_3lCR2 + lnL_Zjets21 + lnL_Zjets22 + lnL_Zjets11 + lnL_Zjets12 + lnL_Zjets01 + lnL_Zjets02 + lnL_WW1 + lnL_WW2;
                                /*cout << i << ")  "
                                << " lnL_all " << lnL_all << " mu_s: " << mu_s << " mu_top " << mu_top << " mu_WZ " << mu_WZ << " mu_Zjets2 " << mu_Zjets2 << " mu_Zjets1 " << mu_Zjets1 << " mu_Zjets0 " << mu_Zjets0 << " mu_WW " << mu_WW << endl;*/
                                if (lnL_all < min_lnL)
                                {
                                    min_lnL = lnL_all;
                                    i_best = i;
                                }
                                // Check if the current iteration count is a multiple of 100,000,000
                                if (i % 1000000 == 0)
                                {
                                    std::cout << "Iteration: " << i << std::endl;
                                }

                                lnL_tot.push_back(lnL_all);
                                mu_s_dummy.push_back(mu_s);
                                mu_top_dummy.push_back(mu_top);
                                mu_WZ_dummy.push_back(mu_WZ);
                                mu_Zjets2_dummy.push_back(mu_Zjets2);
                                mu_Zjets1_dummy.push_back(mu_Zjets1);
                                mu_Zjets0_dummy.push_back(mu_Zjets0);
                                mu_WW_dummy.push_back(mu_WW);
                                i++;
                            }
                        }
                    }
                }
            }
        }
    }
    cout << " the total number of iterations that performed was " << i << endl;
    // Meta tis for loopes
    min = calc_min(lnL_tot, i);

    cout << " --- SIMULTANEOUS FIT: ----" << endl;
    cout << "min is: " << min << "  i_best: " << i_best << " Ln_min: " << lnL_tot[min] << endl;
    cout << "mu_s: " << mu_s_dummy[min] << ", mu_top " << mu_top_dummy[min] << ",  mu_WZ   " << mu_WZ_dummy[min] << ",  muZjets2    " << mu_Zjets2_dummy[min] << ",  muZjets1    " << mu_Zjets1_dummy[min] << ",  muZjets0    " << mu_Zjets0_dummy[min] << ", mu_WW   " << mu_WW_dummy[min] << endl;

    // here is a method to obtain the statistical error for a specific scaling factor each time (need to make it run only for an individual scaling factor)

    // Whichever scaling factor's error i want to find i don't fix it
    mu_WZ = mu_WZ_dummy[min];
    mu_top = mu_top_dummy[min];
    mu_Zjets2 = mu_Zjets2_dummy[min];
    mu_Zjets1 = mu_Zjets1_dummy[min]; // Here needs to be changed depending on what scaling factor you 're working with
    mu_Zjets0 = mu_Zjets0_dummy[min];
    mu_WW = mu_WW_dummy[min];
    // mu_s = mu_s_dummy[min];

    // Plot Likelihood around the minimum:
    i = 0;
    i_best = 0;
    // int n_n = 100;

    for (mu_s = 0.95; mu_s <= 1.05; mu_s += 0.001) // (1.2-0.8 = 0.4/0.001->400) [Here needs to be changed depending on what scaling factor you 're working with]
    {                                              // i want t o calculate the same likelihoods as before so i need this function again calc_neglnL!!
        lnL_signal1 = calc_neglnL(1295, 667.74, mu_s, 78.38, mu_top, 339.55, mu_WZ, 1.22, mu_Zjets2, 39.51, mu_Zjets1, 72.03, mu_Zjets0, 18.52, mu_WW, 21.65);

        lnL_signal2 = calc_neglnL(1426, 939.55, mu_s, 22.94, mu_top, 348.29, mu_WZ, 0.22, mu_Zjets2, 3.14, mu_Zjets1, 35.89, mu_Zjets0, 19.99, mu_WW, 24.32);

        lnL_top1 = calc_neglnL(2867, 0.0, mu_s, 2860.72, mu_top, 0.38, mu_WZ, 0.06, mu_Zjets2, 0.00, mu_Zjets1, 0.0, mu_Zjets0, 6.07, mu_WW, 0.92);

        lnL_top2 = calc_neglnL(2869, 0.0, mu_s, 2792.99, mu_top, 0.48, mu_WZ, 0.10, mu_Zjets2, 0.12, mu_Zjets1, 0.0, mu_Zjets0, 6.96, mu_WW, 3.64);

        lnL_3lCR1 = calc_neglnL(1229, 0.86, mu_s, 34.00, mu_top, 1080.97, mu_WZ, 9.18, mu_Zjets2, 20.52, mu_Zjets1, 2.45, mu_Zjets0, 0.37, mu_WW, 51.30);

        lnL_3lCR2 = calc_neglnL(1180, 1.42, mu_s, 22.97, mu_top, 1018.03, mu_WZ, 6.53, mu_Zjets2, 13.96, mu_Zjets1, 33.29, mu_Zjets0, 0.27, mu_WW, 62.80);

        lnL_Zjets21 = calc_neglnL(2871, 83.64, mu_s, 152.35, mu_top, 158.64, mu_WZ, 2296.70, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 12.44, mu_WW, 24.86);

        lnL_Zjets22 = calc_neglnL(2898, 78.29, mu_s, 91.75, mu_top, 150.05, mu_WZ, 2362.55, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 5.95, mu_WW, 33.03);

        lnL_Zjets11 = calc_neglnL(2836, 86.51, mu_s, 144.00, mu_top, 144.37, mu_WZ, 0.0, mu_Zjets2, 1790.57, mu_Zjets1, 0.0, mu_Zjets0, 34.78, mu_WW, 18.44);

        lnL_Zjets12 = calc_neglnL(2897, 58.53, mu_s, 52.53, mu_top, 103.86, mu_WZ, 0.0, mu_Zjets2, 2007.90, mu_Zjets1, 0.0, mu_Zjets0, 11.20, mu_WW, 27.66);

        lnL_Zjets01 = calc_neglnL(1355, 36.96, mu_s, 63.07, mu_top, 40.47, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 855.74, mu_Zjets0, 19.71, mu_WW, 7.57);

        lnL_Zjets02 = calc_neglnL(1336, 60.77, mu_s, 27.94, mu_top, 38.90, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 874.37, mu_Zjets0, 27.96, mu_WW, 4.84);

        lnL_WW1 = calc_neglnL(1445, -0.08, mu_s, 1067.54, mu_top, 8.99, mu_WZ, 0.19, mu_Zjets2, 2.29, mu_Zjets1, 0.29, mu_Zjets0, 206.19, mu_WW, 28.21);

        lnL_WW2 = calc_neglnL(1458, 0.49, mu_s, 767.31, mu_top, 13.30, mu_WZ, 0.21, mu_Zjets2, -0.07, mu_Zjets1, 0.35, mu_Zjets0, 425.87, mu_WW, 68.67);

        lnL_all = lnL_signal1 + lnL_signal2 + lnL_top1 + lnL_top2 + lnL_3lCR1 + lnL_3lCR2 + lnL_Zjets21 + lnL_Zjets22 + lnL_Zjets11 + lnL_Zjets12 + lnL_Zjets01 + lnL_Zjets02 + lnL_WW1 + lnL_WW2;

        cout << i << ")  "
             << " lnL_all " << lnL_all << " mu_s: " << mu_s << " mu_top " << mu_top << " mu_WZ " << mu_WZ << " mu_Zjets2 " << mu_Zjets2 << " mu_Zjets1 " << mu_Zjets1 << " mu_Zjets0 " << mu_Zjets0 << " mu_WW " << mu_WW << endl;
        if (lnL_all < Min_lnL) // maybe i need to apply <= because the same lnL_all will find as previously
        {
            Min_lnL = lnL_all;
            i_best = i;
        }
        lnL_tot1[i] = lnL_all;
        mu_s_dummy1[i] = mu_s; // Here needs to be changed depending on what scaling factor you 're working with
        i++;
    }

    cout << " --- SIMULTANEOUS FIT: ----" << endl;
    cout << " min is: " << min << " Ln_min: " << lnL_tot[min] << endl;
    cout << " mu_s: " << mu_s_dummy[min] << " mu_top " << mu_top_dummy[min] << ", mu_WZ " << mu_WZ_dummy[min] << ", muZjets2 " << mu_Zjets2_dummy[min] << ", muZjets1 " << mu_Zjets1_dummy[min] << ", muZjets0 " << mu_Zjets0_dummy[min] << ", mu_WW " << mu_WW_dummy[min] << endl;
    mu_WZ = mu_WZ_dummy[min];
    mu_top = mu_top_dummy[min];
    mu_Zjets2 = mu_Zjets2_dummy[min];
    mu_Zjets1 = mu_Zjets1_dummy[min]; // Here needs to be changed depending on what scaling factor you 're working with
    mu_Zjets0 = mu_Zjets0_dummy[min];
    mu_WW = mu_WW_dummy[min];
    // mu_s = mu_s_dummy[min];

    cout << " --- FIT only mu_s , keeing the rest fixed: ----" << endl;
    cout << "min is at i_best: " << i_best << " Ln_min: " << lnL_tot1[i_best] << endl; // Here needs to be changed depending on what scaling factor you 're working with. Underneath
    cout << "mu_s: " << mu_s_dummy1[i_best] << " mu_top " << mu_top << " mu_WZ " << mu_WZ << " mu_Zjets2 " << mu_Zjets2 << " mu_Zjets1 " << mu_Zjets1 << " mu_Zjets0 " << mu_Zjets0 << " mu_WW " << mu_WW << endl;

    double *xValues = new double[i];
    std::copy(mu_s_dummy1.begin(), mu_s_dummy1.end(), xValues); // Here needs to be changed depending on what scaling factor you 're working with
    double *yValues = new double[i];
    std::copy(lnL_tot1.begin(), lnL_tot1.end(), yValues);
    TCanvas *c = new TCanvas("c", "Likelihood", 0., 0., 700, 700);
    TGraph *gr = new TGraph(80, xValues, yValues);
    cout << i_best << endl;
    double dummy;
    for (int q = i_best; q <= i_best + 1000; q++)
    {
        if (lnL_tot1[q] >= lnL_tot1[i_best] + 0.5)
        {
            cout << "error+:" << mu_s_dummy1[q] - mu_s_dummy1[i_best] << endl; // Here needs to be changed depending on what scaling factor you 're working with
            q = i_best + 1000;
        }
    }
    for (int q = i_best; q >= i_best - 1000; q--)
    {
        if (lnL_tot1[q] >= lnL_tot1[i_best] + 0.5)
        {
            cout << "error-:" << mu_s_dummy1[q] - mu_s_dummy1[i_best] << endl; // Here needs to be changed depending on what scaling factor you 're working with
            break;
        }
    }

    // Create a TLine for the horizontal line
    double horizontalLineY = lnL_tot1[i_best] + 0.5;
    TLine *horizontalLine = new TLine(gr->GetXaxis()->GetXmin(), horizontalLineY, gr->GetXaxis()->GetXmax(), horizontalLineY);
    horizontalLine->SetLineColor(kRed);
    horizontalLine->SetLineStyle(2); // You can adjust the line style as needed

    const char *additionalText[] = {
        "#leftarrow error (+) #rightarrow",
        "#leftarrow error (-) #rightarrow"};

    // Create TLines for vertical lines
    double verticalLineX1 = mu_s_dummy1[i_best] - 0.1; // Adjust these values as needed
    TLine *verticalLine1 = new TLine(verticalLineX1, gr->GetYaxis()->GetXmin(), verticalLineX1, gr->GetYaxis()->GetXmax());

    verticalLine1->SetLineColor(kRed);
    verticalLine1->SetLineStyle(2); // You can adjust the line style as needed
    c->cd();
    gr->SetMarkerStyle(10);
    gr->SetTitle("");
    gr->GetXaxis()->SetTitle("#mu_{s}");
    gr->GetYaxis()->SetTitle("-lnL_{tot}");
    gr->Draw();
    horizontalLine->Draw("same");
    verticalLine1->Draw("same");
    // Add lines of text using TLatex
    double linePosition = 0.75; // Initial position for the first line
    for (const char *text : additionalText)
    {
        TLatex *textLabel = new TLatex(0.1, linePosition, text);
        textLabel->SetNDC();
        textLabel->SetTextFont(12);
        textLabel->SetTextSize(0.04);
        textLabel->Draw("same");

        linePosition -= 0.07; // Adjust this value to control the spacing between lines
    }

    // Here is a method to obtain again statistical error on a specific scaling factor (need to make it run only for an individual scaling factor)
    TCanvas *c1 = new TCanvas("c1", "Statistical", 0., 0., 700, 700);
    TH1F *h_mu_s_valid = new TH1F("H. Statistical", "", 160, 0.6, 1.40);
    TRandom3 rnd;

    std::vector<double> LnL_tot2;
    std::vector<double> mu_s_dummy3;

    for (int i = 0; i < 50000; i++)
    {
        int count_valid = 0;
        int Min_valid = 0;
        // int I_best = 0;
        // int count = 0;

        //////////////// mu_s  ////////////////
        Double_t N_SR_DATA1 = rnd.Gaus(1295, sqrt(1295)); // it take values form a Gaus with mean the N_DATA = sum of monte carlo and deviation the sqrt(N_DATA)

        Double_t N_SR_DATA2 = rnd.Gaus(1426, sqrt(1426));

        /////////////// mu_top /////////////////
        Double_t N_emCR_B_DATA1 = rnd.Gaus(2867, sqrt(2867));

        Double_t N_emCR_B_DATA2 = rnd.Gaus(2869, sqrt(2869));

        ////////////// mu_WZ /////////////////
        Double_t N_3lCR_DATA1 = rnd.Gaus(1229, sqrt(1229));

        Double_t N_3lCR_DATA2 = rnd.Gaus(1180, sqrt(1180));

        ////////////// mu_Zjets2 //////////////
        Double_t N_Zjets2_DATA1 = rnd.Gaus(2871, sqrt(2871));

        Double_t N_Zjets2_DATA2 = rnd.Gaus(2898, sqrt(2898));

        ////////////// mu_Zjets1 ///////////////
        Double_t N_Zjets1_DATA1 = rnd.Gaus(2836, sqrt(2836));

        Double_t N_Zjets1_DATA2 = rnd.Gaus(2897, sqrt(2897));

        ////////////// mu_Zjets0 ///////////////
        Double_t N_Zjets0_DATA1 = rnd.Gaus(1355, sqrt(1355));

        Double_t N_Zjets0_DATA2 = rnd.Gaus(1336, sqrt(1336));

        ////////////// mu_WW ///////////////
        Double_t N_emCR_A_DATA1 = rnd.Gaus(1445, sqrt(1445));

        Double_t N_emCR_A_DATA2 = rnd.Gaus(1458, sqrt(1458));

        /*cout << "Outer Iteration " << i << endl;
        cout << "N SR_DATA is : " << ((2721 - N_SR_DATA)/sqrt(2721))  << " sigma " << endl;
        cout << "N emCR_B_DATA is : " << ((5736 - N_emCR_B_DATA)/ 75.74) << " sigma " << endl;*/

        for (mu_s = 0.50; mu_s <= 1.50; mu_s += 0.001) // Here need adjustment (100 values for vectors)
        {
            // cout << "Outer Iteration " << i << " Inner iteration " << count << ", mu_s = " << mu_s << endl;
            // I_best += 1;
            // count += 1;

            lnL_signal1 = calc_neglnL(N_SR_DATA1, 667.74, mu_s, 78.38, mu_top, 339.55, mu_WZ, 1.22, mu_Zjets2, 39.51, mu_Zjets1, 72.03, mu_Zjets0, 18.52, mu_WW, 21.65);

            lnL_signal2 = calc_neglnL(N_SR_DATA2, 939.55, mu_s, 22.94, mu_top, 348.29, mu_WZ, 0.22, mu_Zjets2, 3.14, mu_Zjets1, 35.89, mu_Zjets0, 19.99, mu_WW, 24.32);

            lnL_top1 = calc_neglnL(N_emCR_B_DATA1, 0.0, mu_s, 2860.72, mu_top, 0.38, mu_WZ, 0.06, mu_Zjets2, 0.00, mu_Zjets1, 0.0, mu_Zjets0, 6.07, mu_WW, 0.92);

            lnL_top2 = calc_neglnL(N_emCR_B_DATA2, 0.0, mu_s, 2792.99, mu_top, 0.48, mu_WZ, 0.10, mu_Zjets2, 0.12, mu_Zjets1, 0.0, mu_Zjets0, 6.96, mu_WW, 3.64);

            lnL_3lCR1 = calc_neglnL(N_3lCR_DATA1, 0.86, mu_s, 34.00, mu_top, 1080.97, mu_WZ, 9.18, mu_Zjets2, 20.52, mu_Zjets1, 2.45, mu_Zjets0, 0.37, mu_WW, 51.30);

            lnL_3lCR2 = calc_neglnL(N_3lCR_DATA2, 1.42, mu_s, 22.97, mu_top, 1018.03, mu_WZ, 6.53, mu_Zjets2, 13.96, mu_Zjets1, 33.29, mu_Zjets0, 0.27, mu_WW, 62.80);

            lnL_Zjets21 = calc_neglnL(N_Zjets2_DATA1, 83.64, mu_s, 152.35, mu_top, 158.64, mu_WZ, 2296.70, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 12.44, mu_WW, 24.86);

            lnL_Zjets22 = calc_neglnL(N_Zjets2_DATA2, 78.29, mu_s, 91.75, mu_top, 150.05, mu_WZ, 2362.55, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 5.95, mu_WW, 33.03);

            lnL_Zjets11 = calc_neglnL(N_Zjets1_DATA1, 86.51, mu_s, 144.00, mu_top, 144.37, mu_WZ, 0.0, mu_Zjets2, 1790.57, mu_Zjets1, 0.0, mu_Zjets0, 34.78, mu_WW, 18.44);

            lnL_Zjets12 = calc_neglnL(N_Zjets1_DATA2, 58.53, mu_s, 52.53, mu_top, 103.86, mu_WZ, 0.0, mu_Zjets2, 2007.90, mu_Zjets1, 0.0, mu_Zjets0, 11.20, mu_WW, 27.66);

            lnL_Zjets01 = calc_neglnL(N_Zjets0_DATA1, 36.96, mu_s, 63.07, mu_top, 40.47, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 855.74, mu_Zjets0, 19.71, mu_WW, 7.57);

            lnL_Zjets02 = calc_neglnL(N_Zjets0_DATA2, 60.77, mu_s, 27.94, mu_top, 38.90, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 874.37, mu_Zjets0, 27.96, mu_WW, 4.84);

            lnL_WW1 = calc_neglnL(N_emCR_A_DATA1, -0.08, mu_s, 1067.54, mu_top, 8.99, mu_WZ, 0.19, mu_Zjets2, 2.29, mu_Zjets1, 0.29, mu_Zjets0, 206.19, mu_WW, 28.21);

            lnL_WW2 = calc_neglnL(N_emCR_A_DATA2, 0.49, mu_s, 767.31, mu_top, 13.30, mu_WZ, 0.21, mu_Zjets2, -0.07, mu_Zjets1, 0.35, mu_Zjets0, 425.87, mu_WW, 68.67);

            lnL_all = lnL_signal1 + lnL_signal2 + lnL_top1 + lnL_top2 + lnL_3lCR1 + lnL_3lCR2 + lnL_Zjets21 + lnL_Zjets22 + lnL_Zjets11 + lnL_Zjets12 + lnL_Zjets01 + lnL_Zjets02 + lnL_WW1 + lnL_WW2;

            LnL_tot2.push_back(lnL_all);

            mu_s_dummy3.push_back(mu_s); // Need adjustment
            count_valid += 1;
        }
        c1->cd();
        Min_valid = calc_min(LnL_tot2, count_valid);
        Double_t fill_valid = mu_s_dummy3[Min_valid]; // Here needs adjustment
        h_mu_s_valid->Fill(fill_valid);
        mu_s_dummy3.clear(); // Here needs adjustment
        LnL_tot2.clear();
    }

    h_mu_s_valid->GetXaxis()->SetTitle("#mu_{s}"); // Here needs adjustment
    h_mu_s_valid->GetYaxis()->SetTitle("Entries");
    // h_mu_s->SetFillColor(46); // Zjets0
    // h_mu_s->SetFillColor(7); // Zjets1
    // h_mu_s->SetFillColor(8); // Zjets2
    h_mu_s_valid->SetFillColor(41); // Signal
    // h_mu_s->SetFillColor(38); // Top
    // h_mu_s->SetFillColor(28); // WW
    // h_mu_s->SetFillColor(93); //WZ
    h_mu_s_valid->Draw("hist");
    Double_t systematic_valid = h_mu_s_valid->GetRMS();

    cout << " The statistical ancertainty (in order to validate the stat ancertainties from Min_likelihood +0.5) is: " << systematic_valid << endl;

    ///////////////// Here is a method in order to derive statistical sigmas obtained by above (for the pull diagram for the statistical uncertainty) /////////////////////////////

    // TCanvas *c11 = new TCanvas("c11", "PUll_statistical", 0., 0., 700, 700);
    // TH1F *Pull_stat = new TH1F("Pull Stat.", "", 200, -10, 10);
    // TRandom3 rnd20;
    // // TH1F *h_mu_stat = new TH1F("h_mu_stat", "", 60, 0.75, 1.35);

    // // std::vector<double> statsigma;

    // for (int i = 0; i < 50000; i++)
    // {
    //     Double_t Min_LNL = 9999;
    //     Double_t best = 0;
    //     Double_t k = 0;
    //     Double_t errorplus = 0.0;
    //     Double_t errorminus = 0.0;
    //     bool foundPlus = false;
    //     bool foundMinus = false;
    //     std::vector<double> lnL_tot11;
    //     std::vector<double> mu_s_dummy11;

    //     //////////////// mu_s  ////////////////
    //     Double_t N_SR_DATA1 = rnd20.Gaus(1295, sqrt(1295)); // it take values form a Gaus with mean the N_DATA = sum of monte carlo and deviation the sqrt(N_DATA)

    //     Double_t N_SR_DATA2 = rnd20.Gaus(1426, sqrt(1426));

    //     /////////////// mu_top /////////////////
    //     Double_t N_emCR_B_DATA1 = rnd20.Gaus(2867, sqrt(2867));

    //     Double_t N_emCR_B_DATA2 = rnd20.Gaus(2869, sqrt(2869));

    //     ////////////// mu_WZ /////////////////
    //     Double_t N_3lCR_DATA1 = rnd20.Gaus(1229, sqrt(1229));

    //     Double_t N_3lCR_DATA2 = rnd20.Gaus(1180, sqrt(1180));

    //     ////////////// mu_Zjets2 //////////////
    //     Double_t N_Zjets2_DATA1 = rnd20.Gaus(2871, sqrt(2871));

    //     Double_t N_Zjets2_DATA2 = rnd20.Gaus(2898, sqrt(2898));

    //     ////////////// mu_Zjets1 ///////////////
    //     Double_t N_Zjets1_DATA1 = rnd20.Gaus(2836, sqrt(2836));

    //     Double_t N_Zjets1_DATA2 = rnd20.Gaus(2897, sqrt(2897));

    //     ////////////// mu_Zjets0 ///////////////
    //     Double_t N_Zjets0_DATA1 = rnd20.Gaus(1355, sqrt(1355));

    //     Double_t N_Zjets0_DATA2 = rnd20.Gaus(1336, sqrt(1336));

    //     ////////////// mu_WW ///////////////
    //     Double_t N_emCR_A_DATA1 = rnd20.Gaus(1445, sqrt(1445));

    //     Double_t N_emCR_A_DATA2 = rnd20.Gaus(1458, sqrt(1458));

    //     /*cout << "Outer Iteration " << i << endl;
    //     cout << "N SR_DATA is : " << ((2721 - N_SR_DATA)/sqrt(2721))  << " sigma " << endl;
    //     cout << "N emCR_B_DATA is : " << ((5736 - N_emCR_B_DATA)/ 75.74) << " sigma " << endl;*/

    //     for (mu_s = 0.50; mu_s <= 1.50; mu_s += 0.001) // Here need adjustment (100 values for vectors)
    //     {
    //         // cout << "Outer Iteration " << i << " Inner iteration " << count << ", mu_s = " << mu_s << endl;
    //         // I_best += 1;
    //         // count += 1;

    //         lnL_signal1 = calc_neglnL(N_SR_DATA1, 667.74, mu_s, 78.38, mu_top, 339.55, mu_WZ, 1.22, mu_Zjets2, 39.51, mu_Zjets1, 72.03, mu_Zjets0, 18.52, mu_WW, 21.65);

    //         lnL_signal2 = calc_neglnL(N_SR_DATA2, 939.55, mu_s, 22.94, mu_top, 348.29, mu_WZ, 0.22, mu_Zjets2, 3.14, mu_Zjets1, 35.89, mu_Zjets0, 19.99, mu_WW, 24.32);

    //         lnL_top1 = calc_neglnL(N_emCR_B_DATA1, 0.0, mu_s, 2860.72, mu_top, 0.38, mu_WZ, 0.06, mu_Zjets2, 0.00, mu_Zjets1, 0.0, mu_Zjets0, 6.07, mu_WW, 0.92);

    //         lnL_top2 = calc_neglnL(N_emCR_B_DATA2, 0.0, mu_s, 2792.99, mu_top, 0.48, mu_WZ, 0.10, mu_Zjets2, 0.12, mu_Zjets1, 0.0, mu_Zjets0, 6.96, mu_WW, 3.64);

    //         lnL_3lCR1 = calc_neglnL(N_3lCR_DATA1, 0.86, mu_s, 34.00, mu_top, 1080.97, mu_WZ, 9.18, mu_Zjets2, 20.52, mu_Zjets1, 2.45, mu_Zjets0, 0.37, mu_WW, 51.30);

    //         lnL_3lCR2 = calc_neglnL(N_3lCR_DATA2, 1.42, mu_s, 22.97, mu_top, 1018.03, mu_WZ, 6.53, mu_Zjets2, 13.96, mu_Zjets1, 33.29, mu_Zjets0, 0.27, mu_WW, 62.80);

    //         lnL_Zjets21 = calc_neglnL(N_Zjets2_DATA1, 83.64, mu_s, 152.35, mu_top, 158.64, mu_WZ, 2296.70, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 12.44, mu_WW, 24.86);

    //         lnL_Zjets22 = calc_neglnL(N_Zjets2_DATA2, 78.29, mu_s, 91.75, mu_top, 150.05, mu_WZ, 2362.55, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 5.95, mu_WW, 33.03);

    //         lnL_Zjets11 = calc_neglnL(N_Zjets1_DATA1, 86.51, mu_s, 144.00, mu_top, 144.37, mu_WZ, 0.0, mu_Zjets2, 1790.57, mu_Zjets1, 0.0, mu_Zjets0, 34.78, mu_WW, 18.44);

    //         lnL_Zjets12 = calc_neglnL(N_Zjets1_DATA2, 58.53, mu_s, 52.53, mu_top, 103.86, mu_WZ, 0.0, mu_Zjets2, 2007.90, mu_Zjets1, 0.0, mu_Zjets0, 11.20, mu_WW, 27.66);

    //         lnL_Zjets01 = calc_neglnL(N_Zjets0_DATA1, 36.96, mu_s, 63.07, mu_top, 40.47, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 855.74, mu_Zjets0, 19.71, mu_WW, 7.57);

    //         lnL_Zjets02 = calc_neglnL(N_Zjets0_DATA2, 60.77, mu_s, 27.94, mu_top, 38.90, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 874.37, mu_Zjets0, 27.96, mu_WW, 4.84);

    //         lnL_WW1 = calc_neglnL(N_emCR_A_DATA1, -0.08, mu_s, 1067.54, mu_top, 8.99, mu_WZ, 0.19, mu_Zjets2, 2.29, mu_Zjets1, 0.29, mu_Zjets0, 206.19, mu_WW, 28.21);

    //         lnL_WW2 = calc_neglnL(N_emCR_A_DATA2, 0.49, mu_s, 767.31, mu_top, 13.30, mu_WZ, 0.21, mu_Zjets2, -0.07, mu_Zjets1, 0.35, mu_Zjets0, 425.87, mu_WW, 68.67);

    //         lnL_all = lnL_signal1 + lnL_signal2 + lnL_top1 + lnL_top2 + lnL_3lCR1 + lnL_3lCR2 + lnL_Zjets21 + lnL_Zjets22 + lnL_Zjets11 + lnL_Zjets12 + lnL_Zjets01 + lnL_Zjets02 + lnL_WW1 + lnL_WW2;

    //         lnL_tot11.push_back(lnL_all);


    //         /*if (lnL_all < Min_LNL) // this is the procedure to find minimum instead of calc_min function!!!
    //         {
    //             Min_LNL = lnL_all;
    //             best = k;
    //         }*/
    //         /*lnL_tot11[k] = lnL_all;
    //         mu_s_dummy11[k] = mu_s;*/
    //         // Here needs to be changed depending on what scaling factor you 're working with
    //         lnL_tot11.push_back(lnL_all);
    //         mu_s_dummy11.push_back(mu_s);
    //         k++;
    //     }

    //     best = calc_min(lnL_tot11, k);

    //     for (int q = best; q <= best + 1000; q++)
    //     {
    //         if (lnL_tot11[q] >= lnL_tot11[best] + 0.5)
    //         {
    //             errorplus = mu_s_dummy11[q] - mu_s_dummy11[best];
    //             foundPlus = true;
    //             break;
    //         }
    //     }

    //     for (int q = best; q >= best - 1000; q--)
    //     {
    //         if (lnL_tot11[q] >= lnL_tot11[best] + 0.5)
    //         {
    //             errorminus = mu_s_dummy11[q] - mu_s_dummy11[best];
    //             foundMinus = true;
    //             break;
    //         }
    //     }

    //     // Calculate the mean error only if both values were found
    //     if (foundPlus && foundMinus)
    //     {
    //         Double_t mean_error = (errorplus - errorminus) / 2;
    //         // sstatsigma.push_back(mean_error);
    //         //  cout << "Iteration " << i << " statistical uncertainty: " << statsigma[i] << endl;
    //         Double_t fill_stat = mu_s_dummy11[best];
    //         // cout << "Iteration Number " << i << " statistical mu_s: " << fill_stat << " and the mean error is  " << mean_error << endl;
    //         Pull_stat->Fill((fill_stat - 1.006) /  mean_error);
    //     }
    //     mu_s_dummy11.clear();
    //     lnL_tot11.clear();
    //     // h_mu_stat->Reset();
    // }
    // // Create a Gaussian fit function
    // TF1 *gaussianFitt = new TF1("gaussianFitt", "gaus(0)", -10, 10);
    // gaussianFitt->SetParameter(0, 0); // Set mean to 0
    // gaussianFitt->SetParameter(1, 1); // Set standard deviation to 1

    // c11->cd();
    // Pull_stat->GetYaxis()->SetTitle("Entries");
    // Pull_stat->GetXaxis()->SetTitle("pull_{stat}");
    // // Fit the histogram with the Gaussian function
    // Pull_stat->Fit(gaussianFitt, "R");
    // Pull_stat->Draw("SAME");

    // // Get the fitted mean and standard deviation
    // Double_t fittedMeanstat = gaussianFitt->GetParameter(1);
    // Double_t fittedStdDevstat = gaussianFitt->GetParameter(2);

    // // Print the fitted values
    // std::cout << "Fitted Mean: " << fittedMeanstat << std::endl;
    // std::cout << "Fitted StdDev: " << fittedStdDevstat << std::endl;

    return;
}
