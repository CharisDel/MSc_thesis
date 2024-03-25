#include <TCanvas.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>

void ScalingFactorssimuktaneous()
{
    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "Scaling Factors", 800, 600);

    // Create arrays to store the scaling factors
    double scalingFactors[] = {1.014, 1.009, 1.042, 1.303, 1.355, 1.432, 1.132};
    int numScalingFactors = sizeof(scalingFactors) / sizeof(scalingFactors[0]);

    // Create arrays for positive-right and negative-left errors
    double scalingErrorsRight[] = {0.01, 0.02, 0.02, 0.02, 0.03, 0.085, 0.03};
    double scalingErrorsLeft[] = {0.01, 0.02, 0.02, 0.02, 0.02, 0.085, 0.03};

    // Create a TGraph to draw the scaling factor points
    TGraph *graph = new TGraph(numScalingFactors);

    // Set the x and y coordinates for each point
    for (int i = 0; i < numScalingFactors; ++i)
    {
        graph->SetPoint(i, scalingFactors[i], i + 1); // X-coordinate is scaling factor, Y-coordinate is height
    }

    // Customize the graph appearance (markers, colors, etc.)
    graph->SetMarkerStyle(20);     // Marker style
    graph->SetMarkerSize(1.2);     // Marker size
    graph->SetMarkerColor(kBlack); // Marker color

    const char *additionalText[] = {
        "#mu_{top}",
        "1.01^{ 0.01}_{- 0.01}",
        "#mu_{WZ}",
        "1.01^{ 0.02}_{-0.02}",
        "#mu_{Zjet2}",
        "1.04^{ 0.02}_{- 0.02}",
        "#mu_{Zjet1}",
        "1.30^{ 0.02}_{- 0.02}",
        "#mu_{Zjet0}",
        "1.36^{ 0.03}_{- 0.02}",
        "#mu_{WW}",
        "1.43^{ 0.09}_{- 0.09}",
        "#mu_{signal}",
        "1.13^{ 0.03}_{- 0.03}",
        "#bf{ATLAS} Internal"};

    // Set the axis properties
    graph->GetXaxis()->SetTitle("");
    graph->GetYaxis()->SetTitle("");     // No label on y-axis
    graph->GetYaxis()->SetNdivisions(0); // No divisions on y-axis
    graph->SetName("Scaling Factors");

    TAxis *axis = graph->GetXaxis();
    axis->SetLimits(0., 3.);

    // Draw the graph on the canvas
    graph->Draw("AP");

    // Add horizontal lines for positive-right and negative-left errors
    for (int i = 0; i < numScalingFactors; ++i)
    {
        // Positive-right error line to the right of the point
        TLine *linePos = new TLine(scalingFactors[i], i + 1, scalingFactors[i] + scalingErrorsRight[i], i + 1);
        linePos->SetLineStyle(1); // Dotted line style
        linePos->Draw("same");

        // Negative-left error line to the left of the point
        TLine *lineNeg = new TLine(scalingFactors[i] - scalingErrorsLeft[i], i + 1, scalingFactors[i], i + 1);
        lineNeg->SetLineStyle(1); // Dotted line style
        lineNeg->Draw("same");
    }

    // Add vertical line at x = 1
    TLine *line = new TLine(1, 0, 1, numScalingFactors + 1);
    line->SetLineStyle(2); // Dotted line style
    line->Draw("same");


    double linePosition = 0.85; // Initial position for the first line

    for (const char *text : additionalText)
    {
        TLatex *textLabel = new TLatex(0.1, linePosition, text);
        textLabel->SetNDC();
        textLabel->SetTextFont(12);
        textLabel->SetTextSize(0.04);
        textLabel->Draw();

        linePosition -= 0.05; // Adjust this value to control the spacing between lines
    }

    // Update the canvas
    canvas->Update();
}
