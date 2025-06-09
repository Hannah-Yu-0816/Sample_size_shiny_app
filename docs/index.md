# Sample Size Simulation App

Welcome to the documentation page for the Sample Size Simulation Shiny App.

This tool helps researchers estimate the required sample size for detecting treatment effect heterogeneity under different simulation settings.

------------------------------------------------------------------------

## About the App

-   **Purpose**: Support research design by simulating scenarios where treatment effects vary across individuals.
-   **Target Users**: Clinical trial designers, biostatisticians, epidemiologists.
-   **Model Used**: Causal Forest with Honesty Method.

------------------------------------------------------------------------

## Background and Motivation

In clinical trials, traditional sample size calculations typically target the overall treatment effect (Average Treatment Effect, ATE), assuming that all patients respond similarly to treatment.

However, in reality, treatment effects can vary substantially across individuals.
This heterogeneity in treatment response necessitates a shift from focusing solely on the overall effect to studying heterogeneous treatment effects (HTE).

Despite increasing methodological advances in modeling HTE, a key gap remains: there is limited guidance on how to determine the appropriate sample size required to detect such heterogeneity reliably. Existing sample size formulas are largely designed for average effects and do not account for the increased complexity and statistical demands of identifying subgroups with differential responses.

Our project specifically addresses the following questions: 

- What sample size is required to reliably detect HTE?
- How can simulation-based methods assist in planning studies that account for treatment effect variability?

Through simulation-based approaches, we aim to provide practical guidance for researchers designing studies that seek to uncover treatment effect heterogeneity.

------------------------------------------------------------------------

## How to Use the App Locally

1.  Clone the repository:

    ``` bash
    git clone https://github.com/yourusername/sample-size-simulation.git
    ```

2.  Open the project folder in RStudio.

3.  Install the required R packages.

4.  Open the shinytest.R file.

5.  Launch the app.

------------------------------------------------------------------------

## Output Interpretation

Two indices are defined to show the sample size needed to reveal the heterogeneity: heterogeneity captured by effect modifiers and success ratio of pure effect modifier leaf.

### Heterogeneity Captured by Effect Modifiers

The percentage of heterogeneity captured is calculated as follows:

$$
\frac{\frac{1}{n_1}\sum\limits_{i \in S_1}CATE\left(\prod_{j=1}^{k} X_{ij}=1;\underset{\sim}{X_{-ij}}\right)
\-
\frac{1}{n_0}\sum\limits_{i \in S_0}CATE\left(\prod_{j=1}^{k} X_{ij}=0;\underset{\sim}{X_{-ij}}\right)}
{ground \ truth \ of \ CATE}
\times
100\%
$$

Where: 

- $S_1$ is the subset where all selected effect modifiers $X_{ij} = 1$ - $S_0$ is the subset where all selected effect modifiers $X_{ij} = 0$
- The ground truth of CATE is obtained from the data simulation process.

### Success Ratio of Pure Effect Modifier Leaf

The second metric focuses on tree structure: We measure the success rate of trees that split on the effect modifier among the forest.

For this, we look at each individual tree in the forest and classify it into one of three categories: success, partial success, and unsuccess.

-   Success: The tree successfully isolates the effect modifier/phenotype group.
-   Partial Success: The tree isolates only a subset of the effect modifier/phenotype group.
-   Failure: The tree splits on noise and cannot isolate the true subgroup at all.

We then calculate the proportion of success, partial success, and failure across all trees in the forest.

Here is an example with just two variables: $X_1$ is the effect modifier, and $X_2$ is other risk factor.
The eight trees in the following image represent all possible splitting structures for two binary variables.

<img src="https://github.com/user-attachments/assets/0c0eee6a-470b-4450-a348-e09be1fc0552" alt="Screenshot 2025-04-28 at 11 48 01 PM" width="1149"/>

-   Trees 1 through 4, and Tree 8, are considered success — they isolate effect modifier $X_1$ cleanly
-   Trees 6 and 7 are partial success — they isolate a subset of $X_1$
-   Tree 5 is a failure — it splits only on noise

For more than one effect modifiers, we define success when the tree can find a leaf that seperates all "on" effect modifiers (i.e., equal to 1).
This mimics a phenotype-driven mechanism, where the treatment only works if a certain biological profile is fully expressed.

------------------------------------------------------------------------

## Simulation Process

Here’s how the simulation is structured.

We start by evaluating a grid of sample sizes.
At each sample size, we repeat the simulation 500 times (or the time you choose to run) to ensure results are stable and not due to randomness.

For each simulation run, we follow four steps: 

1. Generate data with effect modifiers we got from the user
2. Build a causal forest using that data
3. Compute the two heterogeneity indices we defined earlier — heterogeneity captured and tree success rate
4. Finally, we record the result, and after 500 runs, we take the median as the final estimate for that sample size

To make this efficient, everything is implemented in parallel.

------------------------------------------------------------------------

## Instructions for Users

### Simulation Controls Panel

The Simulation Controls Panel allows users to define the core structure and settings of the simulation.
These parameters control how the synthetic data is generated and how the causal forest is constructed and applied.
The simulation data generation formula and example data are shown on the right side of the panel.
Below is a step-by-step explanation of each field:

1.  Structure Settings

-   Number of Effect Modifiers: This setting controls how many covariates act as treatment effect modifiers.
    -   The first k covariates (e.g., $X_1$, $X_2$, ..., $X_k$) will be treated as effect modifiers.

        For example:

        -   If set to 2 → $X_1$ and $X_2$ are effect modifiers
        -   If set to 4 → $X_1$, $X_2$, $X_3$ and $X_4$

    -   Maximum allowed: 5

        -   The input will be capped at 5 to limit complexity

    -   Adjusting this value will dynamically update:

        -   The data generation formula displayed on the right panel
        -   The sidebar to show input fields for each effect modifier's parameters
-   Number of Other Risk Factors: Defines how many additional covariates affect baseline risk but not the treatment effect.
    -   These covariates are labeled $X_{k+1}$, $X_{k+2}$, ..., continuing from the last effect modifier.

        For example:

        -   If there are 2 effect modifiers → Other Risk Factors start from $X_3$
        -   If there are 5 effect modifiers → Other Risk Factors start from $X_6$

    -   No upper limit, but increasing this number will add more covariate to the simulation

    -   Updating this value will dynamically:

        -   Extend the data generation formula shown on the right panel
        -   Add new parameter input fields in the sidebar for each added risk factor

2.  Risk and Treatment Parameters

-   Baseline Risk Intercept ($\beta_0$)
    -   Sets the base level of risk when no covariates are present. Default is 0
-   Baseline Treatment Effect Intercept ($\gamma_0$)
    -   Controls the average treatment effect when no heterogeneity is present. Default is 0
-   Treatment Effect Size ($\gamma$)
    -   Determines the stength of heterogeneity introduced by the effect modifiers. Default is 1.0
    -   Larger values -\> greater treatment effect variation across subgroups.

3.  Covariate-Specific Configuration

For each covariate (effect modifier or other risk factor), two input fields will appear in the sidebar:

-   Baseline Covariate Effect ($\beta$)
    -   Impact of $X$ on baseline risk. Default is 0.
-   Probability
    -   Sets the prevalence of each covariate (e.g. 0.5 for a binary variable with equal split). Default is 0.5

4.  Simulation Settings

-   Simulation Iterations
    -   Number of repetitions per parameter setting. More iterations yield more stable results, at the cost of computation time.
    -   Default: 500
-   Number of Trees in the Forest
    -   Controls the number of trees in the causal forest algorithm. More trees usually improve stability and predictive power.
    -   Default: 1000
-   Sample Sizes
    -   Comma-separated list of sample sizes to simulate. The app will loop through each size independently and produce comparative plots.
    -   Example input: 100,200,400,600,800,1000

5.  Run Simulation

After configuring all the inputs above, click the “**Run Simulation**” button to begin.
Progress indicators or plots will appear once computation finishes.
Simulation progress will be shown at the bottom right of the shiny app.

6.  Simulation Results

After clicking "Run Simulation", the application displays simulation results as plots. 
These results help evaluate how well the causal forest captures treatment effect heterogeneity under different sample sizes and simulation settings.

We use the following plot as an example of the results (settings: $X_1$ is an effect modifier, $X_2$ is an other risk factor).

<img width="1021" alt="Screenshot 2025-05-19 at 10 45 38 PM" src="https://github.com/user-attachments/assets/bf529b61-3fb5-4120-a6ee-67da3d4768e1" />

This line chart tracks the performance of the model across varying sample sizes. The x-axis represents the sample size, while the y-axis shows the percentage (%) for each metric. Main metrics included:

- Heterogeneity Captured (black line)
  - Percentage of simulations where the model correctly identified the treatment effect heterogeneity driven by $X_1$
  - Improves with larger sample sizes

- Success Ratio of Pure Effect Modifier Leaf (blue line)
  - Trees that successfully isolate the effect modifier/phenotype
  - Improves with larger sample sizes
 
Users can use this plot to access how sample size impacts HTE and typically we select the "elbow point" or point where the performace curve starts to flatten. 

In the example above:

- Both Heterogeneity Captured and Success X Tree increase rapidly between sample sizes 100 → 400
- Around n = 400, the curves begin to stabilize
- Thus, 400 may be considered the overall sample size for this scenario and will be used as input in the parameter tuning panel

### Parameter Tuning Panel

After selecting the overall sample size for your clinical trials.
Then you can come to the second panel: Parameter Tuning.
In this panel, you will have a chance to tune the sample.fraction and honesty.fraction parameters to see their impact on the two indices.

1. Parameters Tuning

- sample.fraction

  Defines the fraction of the full dataset to use in each tree of the causal forest.

  - Value range: (0, 1]
  - At each tree, only a subset of the data is randomly selected
  - Smaller values increase variance but reduce correlation between trees (more diversity)
  - Important restriction: If honesty method used, then sample.fraction must not exceed 0.5

  Example:

  If your total dataset has 1000 samples and
  - sample.fraction = 0.5 -> random 500 observations are used per tree
  - sample.fraction = 0.8 -> not allowed when honesty method is used

- honesty.fraction

  Controls how that sampled subset is split into training and estimation parts, to ensure honest estimation of treatment effects.

  - Honesty means: one part of the data is used to decide where to split (train), the other part is used to estimate effects (evaluate)
  - Helps reduce overfitting and bias
 
  Example:

  With sample.fraction = 0.5 and honesty.fraction = 0.3:

  - Total samples = 1000
  - → 500 samples used in each tree (due to sample.fraction)
  - → Of those 500 samples:
    - 30% (150 samples) used to train the tree (i.e., determine splits)
    - 70% (350 samples) used to estimate treatment effects at the leaves

2. Run Parameter Tuning Simulation

After setting the overall sample size (from **Simulation Controls Panel**) and parameter ranges, click “**Run Parameter Tunig Simulation**” button to execute the tuning simulation.

What it does:
- Runs multiple simulations across the selected parameter grid
- Compares model performance under different configurations

Progress indicators or plots will appear once computation finishes.

Simulation progress will be shown at the bottom right of the shiny app.

3. Simulation Results

After clicking “Run Parameter Tuning Simulation”, two heatmaps will be generated to help evaluate how different combinations of sample.fraction and honesty.fraction affect model performance.

These results are useful for identifying optimal sampling configurations under a fixed total sample size.

- heatmap: Heterogeneity Captured

  Example:

  The following plot used the same parameters setting in the Simulation Controls Panel and used 400 as the overall sample size. 

  <img width="1023" alt="Screenshot 2025-05-19 at 11 24 46 PM" src="https://github.com/user-attachments/assets/d2a2d330-1224-4d59-a5d2-5647dd0a02f4" />

  - Darker → lower heterogeneity capture
  - Brighter → better model performance

- heatmap: Success Ratio of Pure Effect Modifier Leaf

  Example:

  The following plot used the same parameters setting in the Simulation Controls Panel and used 400 as the overall sample size. 

  <img width="1024" alt="Screenshot 2025-05-19 at 11 25 00 PM" src="https://github.com/user-attachments/assets/f7f33f7f-5ef2-42a3-be2e-5243bd52946c" />

  - Darker → lower heterogeneity capture
  - Brighter → better model performance


