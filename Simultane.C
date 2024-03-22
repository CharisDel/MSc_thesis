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

void Simultane()
{
    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    // different for every area
    long double lnL_signal;
    long double lnL_top;
    long double lnL_3lCR;
    long double lnL_Zjets2;
    long double lnL_Zjets1;
    long double lnL_Zjets0;
    long double lnL_WW;
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

    for (mu_top = 0.99; mu_top <= 1.02; mu_top += 0.002) // boundaries for every scaling factor are set +-5σ from its nominal value that has been found
    {
        for (mu_WZ = 0.99; mu_WZ <= 1.02; mu_WZ += 0.002)
        {
            for (mu_Zjets2 = 1.05; mu_Zjets2 <= 1.08; mu_Zjets2 += 0.002)
            {
                for (mu_Zjets1 = 1.31; mu_Zjets1 <= 1.33; mu_Zjets1 += 0.002)
                {
                    for (mu_Zjets0 = 1.34; mu_Zjets0 <= 1.36; mu_Zjets0 += 0.002)
                    {
                        for (mu_WW = 1.45; mu_WW <= 1.48; mu_WW += 0.002)
                        {
                            for (mu_s = 0.99; mu_s <= 1.01; mu_s += 0.001)
                            {
                                                       // obs  signal         top             WZ             Z2               Z1                Z0                 WW            Other
                                lnL_signal = calc_neglnL(2721, 1607.30, mu_s, 101.32, mu_top, 687.84, mu_WZ, 1.44, mu_Zjets2, 42.65, mu_Zjets1, 107.92, mu_Zjets0, 38.51, mu_WW, 45.97);

                                lnL_top = calc_neglnL(5736, 0.0, mu_s, 5653.70, mu_top, 0.87, mu_WZ, 0.15, mu_Zjets2, 0.12, mu_Zjets1, 0.0, mu_Zjets0, 13.03, mu_WW, 4.56);

                                lnL_3lCR = calc_neglnL(2409, 2.29, mu_s, 56.96, mu_top, 2098.99, mu_WZ, 15.70, mu_Zjets2, 34.48, mu_Zjets1, 35.74, mu_Zjets0, 0.65, mu_WW, 114.10);

                                lnL_Zjets2 = calc_neglnL(5769, 161.94, mu_s, 244.09, mu_top, 308.68, mu_WZ, 4659.26, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 18.39, mu_WW, 57.89);

                                lnL_Zjets1 = calc_neglnL(5733, 145.04, mu_s, 196.52, mu_top, 248.24, mu_WZ, 0.0, mu_Zjets2, 3798.46, mu_Zjets1, 0.0, mu_Zjets0, 45.98, mu_WW, 46.10);

                                lnL_Zjets0 = calc_neglnL(2691, 97.73, mu_s, 91.01, mu_top, 79.37, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 1730.11, mu_Zjets0, 47.68, mu_WW, 12.41);

                                lnL_WW = calc_neglnL(2903, 0.41, mu_s, 1834.84, mu_top, 22.30, mu_WZ, 0.40, mu_Zjets2, 2.22, mu_Zjets1, 0.63, mu_Zjets0, 632.07, mu_WW, 96.88);

                                lnL_all = lnL_signal + lnL_top + lnL_3lCR + lnL_Zjets2 + lnL_Zjets1 + lnL_Zjets0 + lnL_WW;
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
    //mu_s = mu_s_dummy[min];

    // Plot Likelihood around the minimum:
    i = 0;
    i_best = 0;
    // int n_n = 100;

    for (mu_s = 0.95; mu_s <= 1.05; mu_s += 0.001) // (1.2-0.8 = 0.4/0.001->400) [Here needs to be changed depending on what scaling factor you 're working with]
    {                                              // i want t o calculate the same likelihoods as before so i need this function again calc_neglnL!!
        lnL_signal = calc_neglnL(2721, 1607.30, mu_s, 101.32, mu_top, 687.84, mu_WZ, 1.44, mu_Zjets2, 42.65, mu_Zjets1, 107.92, mu_Zjets0, 38.51, mu_WW, 45.97);

        lnL_top = calc_neglnL(5736, 0.0, mu_s, 5653.70, mu_top, 0.87, mu_WZ, 0.15, mu_Zjets2, 0.12, mu_Zjets1, 0.0, mu_Zjets0, 13.03, mu_WW, 4.56);

        lnL_3lCR = calc_neglnL(2409, 2.29, mu_s, 56.96, mu_top, 2098.99, mu_WZ, 15.70, mu_Zjets2, 34.48, mu_Zjets1, 35.74, mu_Zjets0, 0.65, mu_WW, 114.10);

        lnL_Zjets2 = calc_neglnL(5769, 161.94, mu_s, 244.09, mu_top, 308.68, mu_WZ, 4659.26, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 18.39, mu_WW, 57.89);

        lnL_Zjets1 = calc_neglnL(5733, 145.04, mu_s, 196.52, mu_top, 248.24, mu_WZ, 0.0, mu_Zjets2, 3798.46, mu_Zjets1, 0.0, mu_Zjets0, 45.98, mu_WW, 46.10);

        lnL_Zjets0 = calc_neglnL(2691, 97.73, mu_s, 91.01, mu_top, 79.37, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 1730.11, mu_Zjets0, 47.68, mu_WW, 12.41);

        lnL_WW = calc_neglnL(2903, 0.41, mu_s, 1834.84, mu_top, 22.30, mu_WZ, 0.40, mu_Zjets2, 2.22, mu_Zjets1, 0.63, mu_Zjets0, 632.07, mu_WW, 96.88);

        lnL_all = lnL_signal + lnL_top + lnL_3lCR + lnL_Zjets2 + lnL_Zjets1 + lnL_Zjets0 + lnL_WW;

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
    //mu_s = mu_s_dummy[min];

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

    ///////////////// Here is a method in order to derive statistical sigmas obtained by above (for the pull diagram for the statistical uncertainty) /////////////////////////////

    TCanvas *c11 = new TCanvas("c11", "PUll_statistical", 0., 0., 700, 700);
    TH1F *Pull_stat = new TH1F("Pull Stat.", "", 100, -10, 10);
    TRandom3 rnd20;
    // TH1F *h_mu_stat = new TH1F("h_mu_stat", "", 60, 0.75, 1.35);

    // std::vector<double> statsigma;

    for (int i = 0; i < 50000; i++)
    {
        Double_t Min_LNL = 9999;
        Double_t best = 0;
        Double_t k = 0;
        Double_t errorplus = 0.0;
        Double_t errorminus = 0.0;
        bool foundPlus = false;
        bool foundMinus = false;
        std::vector<double> lnL_tot11;
        std::vector<double> mu_s_dummy11;

        //////////////// mu_s  ////////////////
        Double_t N_SR_DATA = rnd20.Gaus(2721, sqrt(2721)); // it take values form a Gaus with mean the N_DATA = sum of monte carlo and deviation the sqrt(N_DATA)

        /////////////// mu_top /////////////////
        Double_t N_emCR_B_DATA = rnd20.Gaus(5736, 75.74);

        ////////////// mu_WZ /////////////////
        Double_t N_3lCR_DATA = rnd20.Gaus(2409, 49.08);

        ////////////// mu_Zjets2 //////////////
        Double_t N_Zjets2_DATA = rnd20.Gaus(5769, 75.95);

        ////////////// mu_Zjets1 ///////////////
        Double_t N_Zjets1_DATA = rnd20.Gaus(5733, 75.72);

        ////////////// mu_Zjets0 ///////////////
        Double_t N_Zjets0_DATA = rnd20.Gaus(2691, 51.87);

        ////////////// mu_WW ///////////////
        Double_t N_emCR_A_DATA = rnd20.Gaus(2903, 53.88);

        for (mu_s = 0.50; mu_s <= 1.50; mu_s += 0.001) // (1.2-0.8 = 0.4/0.001->400) [Here needs to be changed depending on what scaling factor you 're working with]
        {                                              // i want t o calculate the same likelihoods as before so i need this function again calc_neglnL!!
            lnL_signal = calc_neglnL(N_SR_DATA, 1607.30, mu_s, 101.32, mu_top, 687.84, mu_WZ, 1.44, mu_Zjets2, 42.65, mu_Zjets1, 107.92, mu_Zjets0, 38.51, mu_WW, 45.97);

            lnL_top = calc_neglnL(N_emCR_B_DATA, 0.0, mu_s, 5653.70, mu_top, 0.87, mu_WZ, 0.15, mu_Zjets2, 0.12, mu_Zjets1, 0.0, mu_Zjets0, 13.03, mu_WW, 4.56);

            lnL_3lCR = calc_neglnL(N_3lCR_DATA, 2.29, mu_s, 56.96, mu_top, 2098.99, mu_WZ, 15.70, mu_Zjets2, 34.48, mu_Zjets1, 35.74, mu_Zjets0, 0.65, mu_WW, 114.10);

            lnL_Zjets2 = calc_neglnL(N_Zjets2_DATA, 161.94, mu_s, 244.09, mu_top, 308.68, mu_WZ, 4659.26, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 18.39, mu_WW, 57.89);

            lnL_Zjets1 = calc_neglnL(N_Zjets1_DATA, 145.04, mu_s, 196.52, mu_top, 248.24, mu_WZ, 0.0, mu_Zjets2, 3798.46, mu_Zjets1, 0.0, mu_Zjets0, 45.98, mu_WW, 46.10);

            lnL_Zjets0 = calc_neglnL(N_Zjets0_DATA, 97.73, mu_s, 91.01, mu_top, 79.37, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 1730.11, mu_Zjets0, 47.68, mu_WW, 12.41);

            lnL_WW = calc_neglnL(N_emCR_A_DATA, 0.41, mu_s, 1834.84, mu_top, 22.30, mu_WZ, 0.40, mu_Zjets2, 2.22, mu_Zjets1, 0.63, mu_Zjets0, 632.07, mu_WW, 96.88);

            lnL_all = lnL_signal + lnL_top + lnL_3lCR + lnL_Zjets2 + lnL_Zjets1 + lnL_Zjets0 + lnL_WW;

            /*if (lnL_all < Min_LNL) // this is the procedure to find minimum instead of calc_min function!!!
            {
                Min_LNL = lnL_all;
                best = k;
            }*/
            /*lnL_tot11[k] = lnL_all;
            mu_s_dummy11[k] = mu_s;*/
            // Here needs to be changed depending on what scaling factor you 're working with
            lnL_tot11.push_back(lnL_all);
            mu_s_dummy11.push_back(mu_s);
            k++;
        }

        best = calc_min(lnL_tot11, k);

        for (int q = best; q <= best + 1000; q++)
        {
            if (lnL_tot11[q] >= lnL_tot11[best] + 0.5)
            {
                errorplus = mu_s_dummy11[q] - mu_s_dummy11[best];
                foundPlus = true;
                break;
            }
        }

        for (int q = best; q >= best - 1000; q--)
        {
            if (lnL_tot11[q] >= lnL_tot11[best] + 0.5)
            {
                errorminus = mu_s_dummy11[q] - mu_s_dummy11[best];
                foundMinus = true;
                break;
            }
        }

        // Calculate the mean error only if both values were found
        if (foundPlus && foundMinus)
        {
            Double_t mean_error_stat = (errorplus - errorminus) / 2;
            // sstatsigma.push_back(mean_error);
            //  cout << "Iteration " << i << " statistical uncertainty: " << statsigma[i] << endl;
            Double_t fill_stat = mu_s_dummy11[best];
            // cout << "Iteration Number " << i << " statistical mu_s: " << fill_stat << " and the mean error is  " << mean_error << endl;
            Pull_stat->Fill((fill_stat - 1.006) / mean_error_stat);
        }
        mu_s_dummy11.clear();
        lnL_tot11.clear();
        // h_mu_stat->Reset();
    }
    // Create a Gaussian fit function
    TF1 *gaussianFitt = new TF1("gaussianFitt", "gaus(0)", -10, 10);
    gaussianFitt->SetParameter(0, 0); // Set mean to 0
    gaussianFitt->SetParameter(1, 1); // Set standard deviation to 1

    c11->cd();
    Pull_stat->GetYaxis()->SetTitle("Entries");
    Pull_stat->GetXaxis()->SetTitle("pull_{stat}");
    // Fit the histogram with the Gaussian function
    Pull_stat->Fit(gaussianFitt, "R");
    Pull_stat->Draw("SAME");

    // Get the fitted mean and standard deviation
    Double_t fittedMeanstat = gaussianFitt->GetParameter(1);
    Double_t fittedStdDevstat = gaussianFitt->GetParameter(2);

    // Print the fitted values
    std::cout << "Fitted Mean: " << fittedMeanstat << std::endl;
    std::cout << "Fitted StdDev: " << fittedStdDevstat << std::endl;

    // Here is a method to obtain again statistical error on a specific scaling factor (need to make it run only for an individual scaling factor)
    TCanvas *c1 = new TCanvas("c1", "Statistical", 0., 0., 700, 700);
    TH1F *h_mu_s_valid = new TH1F("H. Statistical", "", 80, 0.8, 1.20);
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
        Double_t N_SR_DATA = rnd.Gaus(2721, sqrt(2721)); // it take values form a Gaus with mean the N_DATA = sum of monte carlo and deviation the sqrt(N_DATA)

        /////////////// mu_top /////////////////
        Double_t N_emCR_B_DATA = rnd.Gaus(5736, 75.74);

        ////////////// mu_WZ /////////////////
        Double_t N_3lCR_DATA = rnd.Gaus(2409, 49.08);

        ////////////// mu_Zjets2 //////////////
        Double_t N_Zjets2_DATA = rnd.Gaus(5769, 75.95);

        ////////////// mu_Zjets1 ///////////////
        Double_t N_Zjets1_DATA = rnd.Gaus(5733, 75.72);

        ////////////// mu_Zjets0 ///////////////
        Double_t N_Zjets0_DATA = rnd.Gaus(2691, 51.87);

        ////////////// mu_WW ///////////////
        Double_t N_emCR_A_DATA = rnd.Gaus(2903, 53.88);

        /*cout << "Outer Iteration " << i << endl;
        cout << "N SR_DATA is : " << ((2721 - N_SR_DATA)/sqrt(2721))  << " sigma " << endl;
        cout << "N emCR_B_DATA is : " << ((5736 - N_emCR_B_DATA)/ 75.74) << " sigma " << endl;*/

        for (mu_s = 0.50; mu_s <= 1.50; mu_s += 0.001) // Here need adjustment (100 values for vectors)
        {
            // cout << "Outer Iteration " << i << " Inner iteration " << count << ", mu_s = " << mu_s << endl;
            // I_best += 1;
            // count += 1;

            lnL_signal = calc_neglnL(N_SR_DATA, 1607.30, mu_s, 101.32, mu_top, 687.84, mu_WZ, 1.44, mu_Zjets2, 42.65, mu_Zjets1, 107.92, mu_Zjets0, 38.51, mu_WW, 45.97);

            lnL_top = calc_neglnL(N_emCR_B_DATA, 0.0, mu_s, 5653.70, mu_top, 0.87, mu_WZ, 0.15, mu_Zjets2, 0.12, mu_Zjets1, 0.0, mu_Zjets0, 13.03, mu_WW, 4.56);

            lnL_3lCR = calc_neglnL(N_3lCR_DATA, 2.29, mu_s, 56.96, mu_top, 2098.99, mu_WZ, 15.70, mu_Zjets2, 34.48, mu_Zjets1, 35.74, mu_Zjets0, 0.65, mu_WW, 114.10);

            lnL_Zjets2 = calc_neglnL(N_Zjets2_DATA, 161.94, mu_s, 244.09, mu_top, 308.68, mu_WZ, 4659.26, mu_Zjets2, 0.0, mu_Zjets1, 0.0, mu_Zjets0, 18.39, mu_WW, 57.89);

            lnL_Zjets1 = calc_neglnL(N_Zjets1_DATA, 145.04, mu_s, 196.52, mu_top, 248.24, mu_WZ, 0.0, mu_Zjets2, 3798.46, mu_Zjets1, 0.0, mu_Zjets0, 45.98, mu_WW, 46.10);

            lnL_Zjets0 = calc_neglnL(N_Zjets0_DATA, 97.73, mu_s, 91.01, mu_top, 79.37, mu_WZ, 0.0, mu_Zjets2, 0.0, mu_Zjets1, 1730.11, mu_Zjets0, 47.68, mu_WW, 12.41);

            lnL_WW = calc_neglnL(N_emCR_A_DATA, 0.41, mu_s, 1834.84, mu_top, 22.30, mu_WZ, 0.40, mu_Zjets2, 2.22, mu_Zjets1, 0.63, mu_Zjets0, 632.07, mu_WW, 96.88);

            lnL_all = lnL_signal + lnL_top + lnL_3lCR + lnL_Zjets2 + lnL_Zjets1 + lnL_Zjets0 + lnL_WW;

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

    // Here is a method to obtain the systematic error on a specific scaling factor (need to make it run only for an individual scaling factor)
    TCanvas *c2 = new TCanvas("c2", "Systematic", 0., 0., 700, 700);
    TH1F *h_mu_s = new TH1F("H. Systematic", "", 80, 0.80, 1.20);
    TRandom3 rnd10;

    std::vector<double> mu_s_dummy2;
    std::vector<double> mu_top_dummy2(100);
    std::vector<double> mu_WZ_dummy2(100);
    std::vector<double> mu_WW_dummy2(100);
    std::vector<double> mu_Zjets2_dummy2(100);
    std::vector<double> mu_Zjets1_dummy2(100);
    std::vector<double> mu_Zjets0_dummy2(100);
    std::vector<double> LnL_tot1;

    for (int i = 0; i < 50000; i++)
    {
        int count = 0;
        int Min = 0;
        // std::vector<double> mu_s_dummy2(100);
        // std::vector<double> LnL_tot1(100);
        // int I_best = 0;
        // int count = 0;
        // rnd0.SetSeed(i + 1);

        //////////////// mu_s  ////////////////
        Double_t N_SIGNAL1 = rnd10.Gaus(1607.30, 13.93);
        Double_t N_WZ1 = rnd10.Gaus(687.84, 6.2);
        Double_t N_TOP1 = rnd10.Gaus(101.32, 2.62);
        Double_t N_WW1 = rnd10.Gaus(38.51, 1.09);
        Double_t N_Z01 = rnd10.Gaus(107.92, 24.16);
        Double_t N_Z11 = rnd10.Gaus(42.65, 10.31);
        Double_t N_Z21 = rnd10.Gaus(1.44, 0.76);
        Double_t N_OTHER1 = rnd10.Gaus(45.97, 2.06);

        /////////////// mu_top /////////////////
        Double_t N_SIGNAL2 = rnd10.Gaus(0.0, 0.0);
        Double_t N_WZ2 = rnd10.Gaus(0.87, 0.17);
        Double_t N_TOP2 = rnd10.Gaus(5653.70, 17.61);
        Double_t N_WW2 = rnd10.Gaus(13.03, 0.71);
        Double_t N_Z02 = rnd10.Gaus(0.00, 0.00);
        Double_t N_Z12 = rnd10.Gaus(0.12, 0.05);
        Double_t N_Z22 = rnd10.Gaus(0.15, 0.06);
        Double_t N_OTHER2 = rnd10.Gaus(4.56, 1.98);

        ////////////// mu_WZ /////////////////
        Double_t N_SIGNAL3 = rnd10.Gaus(2.29, 0.54);
        Double_t N_WZ3 = rnd10.Gaus(2098.99, 10.23);
        Double_t N_TOP3 = rnd10.Gaus(56.96, 1.82);
        Double_t N_WW3 = rnd10.Gaus(0.65, 0.14);
        Double_t N_Z03 = rnd10.Gaus(35.74, 8.28);
        Double_t N_Z13 = rnd10.Gaus(34.48, 5.81);
        Double_t N_Z23 = rnd10.Gaus(15.70, 4.09);
        Double_t N_OTHER3 = rnd10.Gaus(114.10, 0.93);

        // rnd0.SetSeed(i + 1);
        ////////////// mu_Zjets2 //////////////
        Double_t N_SIGNAL4 = rnd10.Gaus(161.94, 2.48);
        Double_t N_WZ4 = rnd10.Gaus(308.68, 3.03);
        Double_t N_TOP4 = rnd10.Gaus(244.09, 3.77);
        Double_t N_WW4 = rnd10.Gaus(18.39, 0.79);
        Double_t N_Z04 = rnd10.Gaus(0.0, 0.0);
        Double_t N_Z14 = rnd10.Gaus(0.0, 0.0);
        Double_t N_Z24 = rnd10.Gaus(4659.26, 110.96);
        Double_t N_OTHER4 = rnd10.Gaus(57.89, 2.57);

        ////////////// mu_Zjets1 ///////////////
        Double_t N_SIGNAL5 = rnd10.Gaus(145.04, 4.00);
        Double_t N_WZ5 = rnd10.Gaus(248.24, 3.47);
        Double_t N_TOP5 = rnd10.Gaus(196.52, 3.66);
        Double_t N_WW5 = rnd10.Gaus(45.98, 1.19);
        Double_t N_Z05 = rnd10.Gaus(0.0, 0.0);
        Double_t N_Z15 = rnd10.Gaus(3798.46, 89.32);
        Double_t N_Z25 = rnd10.Gaus(0.0, 0.0);
        Double_t N_OTHER5 = rnd10.Gaus(46.10, 2.52);

        ////////////// mu_Zjets0 ///////////////
        Double_t N_SIGNAL6 = rnd10.Gaus(97.73, 3.59);
        Double_t N_WZ6 = rnd10.Gaus(79.37, 2.01);
        Double_t N_TOP6 = rnd10.Gaus(91.01, 2.68);
        Double_t N_WW6 = rnd10.Gaus(47.68, 1.22);
        Double_t N_Z06 = rnd10.Gaus(1730.11, 99.64);
        Double_t N_Z16 = rnd10.Gaus(0.0, 0.0);
        Double_t N_Z26 = rnd10.Gaus(0.0, 0.0);
        Double_t N_OTHER6 = rnd10.Gaus(12.41, 2.43);

        ////////////// mu_WW ///////////////
        Double_t N_SIGNAL7 = rnd10.Gaus(0.41, 0.20);
        Double_t N_WZ7 = rnd10.Gaus(22.30, 1.02);
        Double_t N_TOP7 = rnd10.Gaus(1834.84, 11.13);
        Double_t N_WW7 = rnd10.Gaus(632.07, 4.55);
        Double_t N_Z07 = rnd10.Gaus(0.63, 0.43);
        Double_t N_Z17 = rnd10.Gaus(2.22, 1.33);
        Double_t N_Z27 = rnd10.Gaus(0.40, 0.36);
        Double_t N_OTHER7 = rnd10.Gaus(96.88, 16.40);

        /*cout << "Outer Iteration " << i << endl;
        cout << "N signal 1 is : " << ((N_SIGNAL1 - 1607.30) / 13.93) << " sigma " << endl;
        cout << "N signal 4 is : " << ((N_SIGNAL4 - 161.94) / 2.48) << " sigma " << endl;
        cout << "N WZ 1 is : " << ((N_WZ1 - 687.84) / 6.2) << " sigma " << endl;
        cout << "N WZ 4 is : " << ((N_WZ4 - 308.68) / 3.03) << " sigma " << endl;*/

        for (mu_s = 0.50; mu_s <= 1.50; mu_s += 0.001) // Here need adjustment (100 values for vectors)
        {
            // cout << "Outer Iteration " << i << " Inner iteration " << count << ", mu_s = " << mu_s << endl;
            // I_best += 1;
            // count += 1;

            lnL_signal = calc_neglnL(2721, N_SIGNAL1, mu_s, N_TOP1, mu_top, N_WZ1, mu_WZ, N_Z21, mu_Zjets2, N_Z11, mu_Zjets1, N_Z01, mu_Zjets0, N_WW1, mu_WW, N_OTHER1);

            lnL_top = calc_neglnL(5736, N_SIGNAL2, mu_s, N_TOP2, mu_top, N_WZ2, mu_WZ, N_Z22, mu_Zjets2, N_Z12, mu_Zjets1, N_Z02, mu_Zjets0, N_WW2, mu_WW, N_OTHER2);

            lnL_3lCR = calc_neglnL(2409, N_SIGNAL3, mu_s, N_TOP3, mu_top, N_WZ3, mu_WZ, N_Z23, mu_Zjets2, N_Z13, mu_Zjets1, N_Z03, mu_Zjets0, N_WW3, mu_WW, N_OTHER3);

            lnL_Zjets2 = calc_neglnL(5769, N_SIGNAL4, mu_s, N_TOP4, mu_top, N_WZ4, mu_WZ, N_Z24, mu_Zjets2, N_Z14, mu_Zjets1, N_Z04, mu_Zjets0, N_WW4, mu_WW, N_OTHER4);

            lnL_Zjets1 = calc_neglnL(5733, N_SIGNAL5, mu_s, N_TOP5, mu_top, N_WZ5, mu_WZ, N_Z25, mu_Zjets2, N_Z15, mu_Zjets1, N_Z05, mu_Zjets0, N_WW5, mu_WW, N_OTHER5);

            lnL_Zjets0 = calc_neglnL(2691, N_SIGNAL6, mu_s, N_TOP6, mu_top, N_WZ6, mu_WZ, N_Z26, mu_Zjets2, N_Z16, mu_Zjets1, N_Z06, mu_Zjets0, N_WW6, mu_WW, N_OTHER6);

            lnL_WW = calc_neglnL(2903, N_SIGNAL7, mu_s, N_TOP7, mu_top, N_WZ7, mu_WZ, N_Z27, mu_Zjets2, N_Z17, mu_Zjets1, N_Z07, mu_Zjets0, N_WW7, mu_WW, N_OTHER7);

            lnL_all = lnL_signal + lnL_top + lnL_3lCR + lnL_Zjets2 + lnL_Zjets1 + lnL_Zjets0 + lnL_WW;

            LnL_tot1.push_back(lnL_all);

            mu_s_dummy2.push_back(mu_s); // Need adjustment
            count += 1;
        }
        c2->cd();
        Min = calc_min(LnL_tot1, count);
        Double_t fill = mu_s_dummy2[Min]; // Here needs adjustment
        h_mu_s->Fill(fill);
        mu_s_dummy2.clear(); // Here needs adjustment
        LnL_tot1.clear();
    }

    h_mu_s->GetXaxis()->SetTitle("#mu_{s}"); // Here needs adjustment
    h_mu_s->GetYaxis()->SetTitle("Entries");
    // h_mu_s->SetFillColor(46); // Zjets0
    // h_mu_s->SetFillColor(7); // Zjets1
    // h_mu_s->SetFillColor(8); // Zjets2
    h_mu_s->SetFillColor(41); // Signal
    // h_mu_s->SetFillColor(38); // Top
    // h_mu_s->SetFillColor(28); // WW
    // h_mu_s->SetFillColor(93); //WZ
    h_mu_s->Draw("hist");
    Double_t systematic = h_mu_s->GetRMS();

    cout << " The systematic ancertainty is: " << systematic << endl;

    /////////////////////////////////Method to vreate pull diagram for the systematic error /////////////////////////////////
    TCanvas *c10 = new TCanvas("c10", "PUll_syst", 0., 0., 700, 700);
    // TH1F *h_mu = new TH1F("h_mu", "", 60, 0.75, 1.35);
    TH1F *Pull_syst = new TH1F("Pull Syst.", "", 100, -10, 10);
    TRandom3 rnd15;

    for (int i = 0; i < 50000; i++)
    {
        int count_syst = 0;
        int Min_syst = 0;
        Double_t errorplusyst = 0.0;
        Double_t errorminusyst = 0.0;
        bool foundPlusyst = false;
        bool foundMinusyst = false;
        std::vector<double> LnL_tot10;
        std::vector<double> mu_s_dummy10;
        // Double_t fill_syst;
        // Double_t systematic_syst;
        // Double_t input;

        //////////////// mu_s  ////////////////
        Double_t N_SIGNAL1 = rnd15.Gaus(1607.30, 13.93);
        Double_t N_WZ1 = rnd15.Gaus(687.84, 6.2);
        Double_t N_TOP1 = rnd15.Gaus(101.32, 2.62);
        Double_t N_WW1 = rnd15.Gaus(38.51, 1.09);
        Double_t N_Z01 = rnd15.Gaus(107.92, 24.16);
        Double_t N_Z11 = rnd15.Gaus(42.65, 10.31);
        Double_t N_Z21 = rnd15.Gaus(1.44, 0.76);
        Double_t N_OTHER1 = rnd15.Gaus(45.97, 2.06);

        /////////////// mu_top /////////////////
        Double_t N_SIGNAL2 = rnd15.Gaus(0.0, 0.0);
        Double_t N_WZ2 = rnd15.Gaus(0.87, 0.17);
        Double_t N_TOP2 = rnd15.Gaus(5653.70, 17.61);
        Double_t N_WW2 = rnd15.Gaus(13.03, 0.71);
        Double_t N_Z02 = rnd15.Gaus(0.00, 0.00);
        Double_t N_Z12 = rnd15.Gaus(0.12, 0.05);
        Double_t N_Z22 = rnd15.Gaus(0.15, 0.06);
        Double_t N_OTHER2 = rnd15.Gaus(4.56, 1.98);

        ////////////// mu_WZ /////////////////
        Double_t N_SIGNAL3 = rnd15.Gaus(2.29, 0.54);
        Double_t N_WZ3 = rnd15.Gaus(2098.99, 10.23);
        Double_t N_TOP3 = rnd15.Gaus(56.96, 1.82);
        Double_t N_WW3 = rnd15.Gaus(0.65, 0.14);
        Double_t N_Z03 = rnd15.Gaus(35.74, 8.28);
        Double_t N_Z13 = rnd15.Gaus(34.48, 5.81);
        Double_t N_Z23 = rnd15.Gaus(15.70, 4.09);
        Double_t N_OTHER3 = rnd15.Gaus(114.10, 0.93);

        // rnd0.SetSeed(i + 1);
        ////////////// mu_Zjets2 //////////////
        Double_t N_SIGNAL4 = rnd15.Gaus(161.94, 2.48);
        Double_t N_WZ4 = rnd15.Gaus(308.68, 3.03);
        Double_t N_TOP4 = rnd15.Gaus(244.09, 3.77);
        Double_t N_WW4 = rnd15.Gaus(18.39, 0.79);
        Double_t N_Z04 = rnd15.Gaus(0.0, 0.0);
        Double_t N_Z14 = rnd15.Gaus(0.0, 0.0);
        Double_t N_Z24 = rnd15.Gaus(4659.26, 110.96);
        Double_t N_OTHER4 = rnd15.Gaus(57.89, 2.57);

        ////////////// mu_Zjets1 ///////////////
        Double_t N_SIGNAL5 = rnd15.Gaus(145.04, 4.00);
        Double_t N_WZ5 = rnd15.Gaus(248.24, 3.47);
        Double_t N_TOP5 = rnd15.Gaus(196.52, 3.66);
        Double_t N_WW5 = rnd15.Gaus(45.98, 1.19);
        Double_t N_Z05 = rnd15.Gaus(0.0, 0.0);
        Double_t N_Z15 = rnd15.Gaus(3798.46, 89.32);
        Double_t N_Z25 = rnd15.Gaus(0.0, 0.0);
        Double_t N_OTHER5 = rnd15.Gaus(46.10, 2.52);

        ////////////// mu_Zjets0 ///////////////
        Double_t N_SIGNAL6 = rnd15.Gaus(97.73, 3.59);
        Double_t N_WZ6 = rnd15.Gaus(79.37, 2.01);
        Double_t N_TOP6 = rnd15.Gaus(91.01, 2.68);
        Double_t N_WW6 = rnd15.Gaus(47.68, 1.22);
        Double_t N_Z06 = rnd15.Gaus(1730.11, 99.64);
        Double_t N_Z16 = rnd15.Gaus(0.0, 0.0);
        Double_t N_Z26 = rnd15.Gaus(0.0, 0.0);
        Double_t N_OTHER6 = rnd15.Gaus(12.41, 2.43);

        ////////////// mu_WW ///////////////
        Double_t N_SIGNAL7 = rnd15.Gaus(0.41, 0.20);
        Double_t N_WZ7 = rnd15.Gaus(22.30, 1.02);
        Double_t N_TOP7 = rnd15.Gaus(1834.84, 11.13);
        Double_t N_WW7 = rnd15.Gaus(632.07, 4.55);
        Double_t N_Z07 = rnd15.Gaus(0.63, 0.43);
        Double_t N_Z17 = rnd15.Gaus(2.22, 1.33);
        Double_t N_Z27 = rnd15.Gaus(0.40, 0.36);
        Double_t N_OTHER7 = rnd15.Gaus(96.88, 16.40);

        /*cout << "Outer Iteration " << i << endl;
        cout << "N signal 1 is : " << ((N_SIGNAL1 - 1607.30) / 13.93) << " sigma " << endl;
        cout << "N signal 4 is : " << ((N_SIGNAL4 - 161.94) / 2.48) << " sigma " << endl;
        cout << "N WZ 1 is : " << ((N_WZ1 - 687.84) / 6.2) << " sigma " << endl;
        cout << "N WZ 4 is : " << ((N_WZ4 - 308.68) / 3.03) << " sigma " << endl;
        cout << " value of Nsignal in 3lcr: " << N_SIGNAL3 << endl;*/

        for (mu_s = 0.50; mu_s <= 1.50; mu_s += 0.001) // Here need adjustment (100 values for vectors)
        {
            // cout << "Outer Iteration " << i << " Inner iteration " << count << ", mu_s = " << mu_s << endl;
            // I_best += 1;
            // count += 1;

            lnL_signal = calc_neglnL(2721, N_SIGNAL1, mu_s, N_TOP1, mu_top, N_WZ1, mu_WZ, N_Z21, mu_Zjets2, N_Z11, mu_Zjets1, N_Z01, mu_Zjets0, N_WW1, mu_WW, N_OTHER1);

            lnL_top = calc_neglnL(5736, N_SIGNAL2, mu_s, N_TOP2, mu_top, N_WZ2, mu_WZ, N_Z22, mu_Zjets2, N_Z12, mu_Zjets1, N_Z02, mu_Zjets0, N_WW2, mu_WW, N_OTHER2);

            lnL_3lCR = calc_neglnL(2409, N_SIGNAL3, mu_s, N_TOP3, mu_top, N_WZ3, mu_WZ, N_Z23, mu_Zjets2, N_Z13, mu_Zjets1, N_Z03, mu_Zjets0, N_WW3, mu_WW, N_OTHER3);

            lnL_Zjets2 = calc_neglnL(5769, N_SIGNAL4, mu_s, N_TOP4, mu_top, N_WZ4, mu_WZ, N_Z24, mu_Zjets2, N_Z14, mu_Zjets1, N_Z04, mu_Zjets0, N_WW4, mu_WW, N_OTHER4);

            lnL_Zjets1 = calc_neglnL(5733, N_SIGNAL5, mu_s, N_TOP5, mu_top, N_WZ5, mu_WZ, N_Z25, mu_Zjets2, N_Z15, mu_Zjets1, N_Z05, mu_Zjets0, N_WW5, mu_WW, N_OTHER5);

            lnL_Zjets0 = calc_neglnL(2691, N_SIGNAL6, mu_s, N_TOP6, mu_top, N_WZ6, mu_WZ, N_Z26, mu_Zjets2, N_Z16, mu_Zjets1, N_Z06, mu_Zjets0, N_WW6, mu_WW, N_OTHER6);

            lnL_WW = calc_neglnL(2903, N_SIGNAL7, mu_s, N_TOP7, mu_top, N_WZ7, mu_WZ, N_Z27, mu_Zjets2, N_Z17, mu_Zjets1, N_Z07, mu_Zjets0, N_WW7, mu_WW, N_OTHER7);

            lnL_all = lnL_signal + lnL_top + lnL_3lCR + lnL_Zjets2 + lnL_Zjets1 + lnL_Zjets0 + lnL_WW;

            LnL_tot10.push_back(lnL_all);

            mu_s_dummy10.push_back(mu_s); // Need adjustment
            count_syst += 1;
        }
        Min_syst = calc_min(LnL_tot10, count_syst);

        for (int q = Min_syst; q <= Min_syst + 1000; q++)
        {
            if (LnL_tot10[q] >= LnL_tot10[Min_syst] + 0.5)
            {
                errorplusyst = mu_s_dummy10[q] - mu_s_dummy10[Min_syst];
                foundPlusyst = true;
                break;
            }
        }

        for (int q = Min_syst; q >= Min_syst - 1000; q--)
        {
            if (LnL_tot10[q] >= LnL_tot10[Min_syst] + 0.5)
            {
                errorminusyst = mu_s_dummy10[q] - mu_s_dummy10[Min_syst];
                foundMinusyst = true;
                break;
            }
        }

        // Calculate the mean error only if both values were found
        if (foundPlusyst && foundMinusyst)
        {
            Double_t mean_error_syst = (errorplusyst - errorminusyst) / 2;
            // sstatsigma.push_back(mean_error);
            //  cout << "Iteration " << i << " statistical uncertainty: " << statsigma[i] << endl;
            Double_t fill_syst = mu_s_dummy10[Min_syst];
            // cout << "Iteration Number " << i << " statistical mu_s: " << fill_stat << " and the mean error is  " << mean_error << endl;
            Pull_syst->Fill((fill_syst - 1.006) / (0.76 * mean_error_syst));
            if(i % 5000 == 0){
                cout << "mean error please be higher than 0.24         " << mean_error_syst << endl;
            }
        }
        mu_s_dummy10.clear();
        LnL_tot10.clear();
    }
    // Create a Gaussian fit function
    TF1 *gaussianFit = new TF1("gaussianFit", "gaus", -10, 10);
    gaussianFit->SetParameter(0, 0); // Set mean to 0
    gaussianFit->SetParameter(1, 1); // Set standard deviation to 1

    c10->cd();

    Pull_syst->GetYaxis()->SetTitle("Entries");
    Pull_syst->GetXaxis()->SetTitle("pull_{syst}");
    // Fit the histogram with the Gaussian function
    Pull_syst->Fit(gaussianFit, "R");
    Pull_syst->Draw("SAME");
    // Get the fitted mean and standard deviation
    Double_t fittedMean = gaussianFit->GetParameter(1);
    Double_t fittedStdDev = gaussianFit->GetParameter(2);

    // Print the fitted values
    std::cout << "Fitted Mean: " << fittedMean << std::endl;
    std::cout << "Fitted StdDev: " << fittedStdDev << std::endl;

    auto end1 = std::chrono::high_resolution_clock::now();
    // Calculate the duration
    auto duration1 = std::chrono::duration_cast<std::chrono::microseconds>(end1 - start).count();

    std::cout << "Time taken by script to run: " << duration1 / pow(10, 6) << " seconds"
              << " or " << duration1 / (60 * pow(10, 6)) << " minutes" << std::endl;

    return;
}
