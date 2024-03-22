#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>

void ScalingFactorsDisplay()
{
    // Create a canvas
    TCanvas *canvas = new TCanvas("canvas", "Scaling Factors", 800, 600);

    // Create arrays to store the scaling factors and their errors
    double scalingFactors[] = {1.012, 1.438, 1.012, 1.066, 1.322, 1.352, 1.006};
    double scalingFactorErrors[] = {0.012, 0.087, 0.022, 0.017, 0.020, 0.03, 0.032};
    int numScalingFactors = sizeof(scalingFactors) / sizeof(scalingFactors[0]);

    // Create a TGraphErrors to draw the scaling factor points with errors
    TGraphErrors *graph = new TGraphErrors(numScalingFactors);

    // Set the x and y coordinates for each point
    for (int i = 0; i < numScalingFactors; ++i)
    {
        graph->SetPoint(i, scalingFactors[i], i + 1);       // X-coordinate is scaling factor, Y-coordinate is height
        graph->SetPointError(i, scalingFactorErrors[i], 0); // Horizontal error is scaling factor error, vertical error is zero
    }

    // Customize the graph appearance (markers, colors, etc.)
    graph->SetMarkerStyle(20);     // Marker style
    graph->SetMarkerSize(1.2);     // Marker size
    graph->SetMarkerColor(kBlack); // Marker color

    const char *additionalText[] = {
        "#mu_{top}",
        "1.012^{ 0.012}_{- 0.012}",
        "#mu_{WW}",
        "1.438^{ 0.087}_{- 0.087}",
        "#mu_{WZ}",
        "1.012^{ 0.022}_{-0.022}",
        "#mu_{Zjet2}",
        "1.066^{ 0.017}_{- 0.017}",
        "#mu_{Zjet1}",
        "1.322^{ 0.020}_{- 0.020}",
        "#mu_{Zjet0}",
        "1.352^{ 0.030}_{- 0.030}",
        "#mu_{signal}",
        "1.006^{ 0.032}_{- 0.032}"};

    // Set the axis properties
    graph->GetXaxis()->SetTitle("");
    graph->GetYaxis()->SetTitle("");     // No label on y-axis
    graph->GetYaxis()->SetNdivisions(0); // No divisions on y-axis
    graph->SetName("Scaling Factors");

    TAxis *axis = graph->GetXaxis();

    axis->SetLimits(0., 3.);         

    // Draw the graph on the canvas
    graph->Draw("AP");

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
