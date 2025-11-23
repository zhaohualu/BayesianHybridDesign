# UI.R

library(shiny)
library(ggplot2)
library(foreach)
library(doParallel)
library(SAMprior)
library(shinyjs)
library(BayesianHybridDesign)
library(future)
library(promises)
library(plotly)
library(DT)
library(shinydashboard)
library(shinyWidgets)

shinyUI(
  dashboardPage(
    skin = "blue",

    # Header
    dashboardHeader(
      title = tags$div(
        tags$img(src = "https://i.ibb.co/yFzV7v9r/IDSWG-logo.png", height = "50px", style = "margin-right: 10px;"),
        "Bayesian Hybrid Design"
      ),
      titleWidth = 350
    ),

    # Sidebar
    dashboardSidebar(
      width = 250,
      sidebarMenu(
        id = "sidebar",
        menuItem("Home", tabName = "home", icon = icon("home")),
        menuItem("Dynamic Power Prior (DPP)", icon = icon("chart-line"),
                 menuSubItem("Library", tabName = "dpp_library", icon = icon("book")),
                 menuSubItem("Single Design", tabName = "dpp_design", icon = icon("calculator")),
                 menuSubItem("Comparative Analysis", tabName = "dpp_table", icon = icon("table")),
                 menuSubItem("Statistical Analysis", tabName = "dpp_analysis", icon = icon("microscope"))
        ),
        menuItem("SAM Prior", icon = icon("layer-group"),
                 menuSubItem("Library", tabName = "sam_library", icon = icon("book")),
                 menuSubItem("Single Design", tabName = "sam_design", icon = icon("calculator")),
                 menuSubItem("Comparative Analysis", tabName = "sam_table", icon = icon("table"))
        ),
        menuItem("Fisher's Exact", icon = icon("dice"),
                 menuSubItem("Library", tabName = "fisher_library", icon = icon("book")),
                 menuSubItem("Power Analysis", tabName = "fisher_power", icon = icon("bolt")),
                 menuSubItem("Bound Analysis", tabName = "fisher_bound", icon = icon("border-all"))
        )
      )
    ),

    # Body
    dashboardBody(
      useShinyjs(),

      tags$head(
        tags$style(HTML("
          /* Main Layout */
          .content-wrapper { background-color: #f8f9fc; }
          .main-header .logo { background-color: #667eea !important; }
          .main-header .navbar { background-color: #667eea !important; }
          .main-sidebar { background-color: #2c3e50 !important; }

          /* Sidebar Styling */
          .sidebar-menu > li > a { color: #ecf0f1 !important; }
          .sidebar-menu > li:hover > a, .sidebar-menu > li.active > a {
            background-color: #34495e !important;
            border-left: 4px solid #667eea;
          }
          .sidebar-menu .treeview-menu > li > a { color: #bdc3c7 !important; }
          .sidebar-menu .treeview-menu > li:hover > a { color: white !important; }

          /* Box Styling */
          .box {
            box-shadow: 0 4px 15px rgba(0,0,0,0.08);
            border-radius: 8px;
            border-top: none;
            margin-bottom: 25px;
          }
          .box-header {
            background: linear-gradient(135deg, #f8f9fc 0%, #ffffff 100%);
            border-bottom: 2px solid #e8e9ef;
            border-radius: 8px 8px 0 0;
            padding: 20px;
          }
          .box-title { font-weight: 600; color: #667eea; font-size: 1.3em; }
          .box-body { padding: 25px; background-color: white; }

          /* Info/Value Boxes */
          .info-box {
            box-shadow: 0 4px 15px rgba(0,0,0,0.08);
            border-radius: 8px;
            border-left: 4px solid #667eea;
          }
          .small-box {
            border-radius: 8px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.08);
          }
          .small-box h3 { font-size: 2.5em; font-weight: 700; }

          /* Buttons */
          .btn-primary {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            border: none;
            border-radius: 25px;
            padding: 12px 35px;
            font-weight: 500;
            transition: all 0.3s;
            box-shadow: 0 4px 10px rgba(102, 126, 234, 0.3);
          }
          .btn-primary:hover {
            transform: translateY(-2px);
            box-shadow: 0 6px 15px rgba(102, 126, 234, 0.4);
          }
          .btn-success {
            background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
            border: none;
            border-radius: 25px;
            padding: 12px 35px;
            font-weight: 500;
            box-shadow: 0 4px 10px rgba(17, 153, 142, 0.3);
          }
          .btn-warning {
            background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
            border: none;
            border-radius: 25px;
            padding: 12px 35px;
            font-weight: 500;
            box-shadow: 0 4px 10px rgba(240, 147, 251, 0.3);
          }

          /* Form Controls */
          .form-control, .selectize-input {
            border-radius: 6px;
            border: 1.5px solid #e0e0e0;
            padding: 10px 15px;
            transition: border 0.3s;
          }
          .form-control:focus, .selectize-input.focus {
            border-color: #667eea;
            box-shadow: 0 0 0 0.2rem rgba(102, 126, 234, 0.15);
          }

          /* Nav Tabs */
          .nav-tabs-custom {
            box-shadow: 0 4px 15px rgba(0,0,0,0.08);
            border-radius: 8px;
            margin-bottom: 25px;
          }
          .nav-tabs-custom > .nav-tabs > li.active > a {
            border-top: 3px solid #667eea;
            font-weight: 600;
            color: #667eea;
          }
          .nav-tabs-custom > .nav-tabs > li > a {
            border-radius: 8px 8px 0 0;
            transition: all 0.3s;
          }
          .nav-tabs-custom > .nav-tabs > li > a:hover {
            background-color: #f8f9fc;
          }

          /* Tables */
          .dataTables_wrapper .dataTables_length select,
          .dataTables_wrapper .dataTables_filter input {
            border-radius: 6px;
            border: 1.5px solid #e0e0e0;
          }
          table.dataTable thead th {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            font-weight: 600;
            border: none;
          }
          table.dataTable tbody tr:hover {
            background-color: #f8f9fc;
          }

          /* Well Panels */
          .well {
            background-color: #f8f9fc;
            border: 1.5px solid #e8e9ef;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
          }

          /* Status Messages */
          .status-message {
            padding: 15px;
            margin: 15px 0;
            border-radius: 8px;
            font-size: 1.05em;
          }

          /* Content Headers */
          .content-header h1 {
            color: #667eea;
            font-weight: 600;
            font-size: 2.2em;
          }

          /* Definition Box */
          .definition-box {
            background-color: #fff;
            padding: 25px;
            border-radius: 8px;
            border-left: 5px solid #667eea;
            box-shadow: 0 2px 8px rgba(0,0,0,0.04);
            min-height: 150px;
          }

          /* Term Links */
          .term-link {
            color: #667eea;
            cursor: pointer;
            padding: 10px 15px;
            display: block;
            border-radius: 6px;
            transition: all 0.3s;
            margin: 5px 0;
          }
          .term-link:hover {
            background-color: #f8f9fc;
            color: #764ba2;
            transform: translateX(5px);
          }

          /* Value Boxes */
          .small-box.bg-blue { background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important; }
          .small-box.bg-yellow { background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%) !important; }
          .small-box.bg-green { background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%) !important; }
          .small-box.bg-orange { background: linear-gradient(135deg, #ff9a44 0%, #fc6076 100%) !important; }

          /* Page Title Styling */
          .content-header {
            margin-bottom: 30px;
            padding: 20px 15px;
            background: linear-gradient(135deg, #f8f9fc 0%, #ffffff 100%);
            border-radius: 8px;
          }

          /* Help Text */
          .help-block {
            color: #7f8c8d;
            font-size: 0.95em;
            margin-top: 8px;
          }

          /* Collapsible Box Animation */
          .box.collapsed-box { transition: all 0.3s; }

          /* Plotly Graphs */
          .plotly { border-radius: 8px; }

          /* Loading Spinner */
          .shiny-output-error { color: #e74c3c; }
          .shiny-output-error:before { content: '⚠ '; }
        "))
      ),

      tabItems(
        # Home Tab
        tabItem(
          tabName = "home",

          # Hero Section
          fluidRow(
            column(12,
                   div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                        padding: 80px 20px;
                        border-radius: 10px;
                        margin-bottom: 40px;
                        box-shadow: 0 10px 30px rgba(0,0,0,0.2);",
                       div(style = "text-align: center; color: white;",
                           tags$img(src = "https://i.postimg.cc/59rbFS2r/design1.png",
                                    height = "140px",
                                    style = "filter: drop-shadow(0 4px 6px rgba(0,0,0,0.3)); margin-bottom: 25px;"),
                           h1(strong("Bayesian Hybrid Design"),
                              style = "font-size: 3.5em; margin-bottom: 20px; text-shadow: 2px 2px 4px rgba(0,0,0,0.3);"),
                           h3("An R Shiny App for Bayesian Hybrid Design and Analysis",
                              style = "font-weight: 300; margin-bottom: 15px; opacity: 0.95; font-size: 1.5em;"),
                           p("Accelerating drug development through intelligent data integration",
                             style = "font-size: 1.1em; opacity: 0.9; margin-top: 20px;")
                       )
                   )
            )
          ),

          # Introduction Section
          fluidRow(
            column(12,
                   div(style = "background: white; padding: 40px; border-radius: 10px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 40px;",
                       h2(icon("info-circle"), " Introduction",
                          style = "color: #667eea; margin-bottom: 25px; font-weight: 600; font-size: 2em;"),
                       p(style = "font-size: 1.15em; color: #555; line-height: 1.9; margin-bottom: 20px;",
                         "In the pharmaceutical industry, developing a new drug is a high-cost, long-term process. For example, the average cost to develop a cancer drug is approximately ",
                         strong("$1.2 billion"), ", and the timeline from initial laboratory research to patient use can exceed ",
                         strong("10 years"), "."),
                       p(style = "font-size: 1.15em; color: #555; line-height: 1.9; margin-bottom: 20px;",
                         "To accelerate this process, one promising approach is the ",
                         strong("hybrid design method"),
                         ". This method allows for the appropriate incorporation of external data into a current study, which enhances go/no-go decision-making for the next phases of development."),
                       p(style = "font-size: 1.15em; color: #555; line-height: 1.9; margin-bottom: 20px;",
                         "Multiple methods have been proposed within Bayesian hybrid designs. A particularly important method, recently proposed by Lu et al (2025), has demonstrated both strong performance and fast computation. However, without a visual graphical interface, it remains challenging for users to apply this method practically."),
                       div(style = "background: linear-gradient(135deg, #e8f4f8 0%, #f0f8ff 100%);
                            padding: 25px;
                            border-radius: 8px;
                            border-left: 5px solid #667eea;
                            margin-top: 25px;",
                           p(style = "font-size: 1.15em; color: #2c3e50; line-height: 1.9; margin: 0;",
                             icon("check-circle", style = "color: #11998e; margin-right: 8px;"),
                             strong("In response, "),
                             "we developed this R Shiny app to provide a user-friendly GUI. It implements the ",
                             strong("Dynamic Power Prior (DPP)"), " and ", strong("SAM Prior"),
                             " methods in Bayesian hybrid design. For benchmarking, we have also included the ",
                             strong("Fisher's Exact"), " approach.")
                       )
                   )
            )
          ),

          # Key Concepts Section
          fluidRow(
            column(12,
                   div(style = "text-align: center; margin: 50px 0 30px 0;",
                       h2("Key Concepts",
                          style = "color: #667eea; margin-bottom: 15px; font-weight: 600; font-size: 2.2em;"),
                       p("Understanding the methods available in this application",
                         style = "font-size: 1.2em; color: #666;")
                   )
            )
          ),

          # Feature Cards
          fluidRow(
            column(4,
                   div(style = "background: white;
                        padding: 40px 30px;
                        border-radius: 10px;
                        box-shadow: 0 4px 15px rgba(0,0,0,0.1);
                        transition: transform 0.3s, box-shadow 0.3s;
                        margin: 10px;
                        min-height: 400px;
                        border-top: 5px solid #667eea;",
                       div(style = "text-align: center; margin-bottom: 25px;",
                           icon("database", style = "font-size: 3.5em; color: #667eea;")
                       ),
                       h3("Dynamic Power Prior (DPP)",
                          style = "color: #667eea; text-align: center; margin-bottom: 20px; font-size: 1.5em;"),
                       p(style = "color: #666; text-align: left; line-height: 1.8; font-size: 1.05em;",
                         "A method that allows for the incorporation of historical data into a current study. It dynamically adjusts the weight given to the historical data based on the similarity between the historical and current control groups. This prevents potential bias from conflicting data while maximizing statistical power.")
                   )
            ),
            column(4,
                   div(style = "background: white;
                        padding: 40px 30px;
                        border-radius: 10px;
                        box-shadow: 0 4px 15px rgba(0,0,0,0.1);
                        transition: transform 0.3s, box-shadow 0.3s;
                        margin: 10px;
                        min-height: 400px;
                        border-top: 5px solid #764ba2;",
                       div(style = "text-align: center; margin-bottom: 25px;",
                           icon("layer-group", style = "font-size: 3.5em; color: #764ba2;")
                       ),
                       h3("Self-Adapting Mixture (SAM) Prior",
                          style = "color: #764ba2; text-align: center; margin-bottom: 20px; font-size: 1.5em;"),
                       p(style = "color: #666; text-align: left; line-height: 1.8; font-size: 1.05em;",
                         "This method creates a prior distribution by mixing an informative prior (based on historical data) with a non-informative, or skeptical, prior. The 'self-adapting' feature means it adjusts the borrowing weight dynamically based on conflicts between new and historical data, providing robust protection against prior-data conflicts.")
                   )
            ),
            column(4,
                   div(style = "background: white;
                        padding: 40px 30px;
                        border-radius: 10px;
                        box-shadow: 0 4px 15px rgba(0,0,0,0.1);
                        transition: transform 0.3s, box-shadow 0.3s;
                        margin: 10px;
                        min-height: 400px;
                        border-top: 5px solid #f093fb;",
                       div(style = "text-align: center; margin-bottom: 25px;",
                           icon("dice", style = "font-size: 3.5em; color: #f093fb;")
                       ),
                       h3("Fisher's Exact Method",
                          style = "color: #f093fb; text-align: center; margin-bottom: 20px; font-size: 1.5em;"),
                       p(style = "color: #666; text-align: left; line-height: 1.8; font-size: 1.05em;",
                         "A statistical test used to analyze the association between two categorical variables, such as treatment response. It is often used as a benchmark for comparison in study design. This non-parametric test provides exact p-values and is particularly useful for small sample sizes.")
                   )
            )
          ),

          # How It Works Section
          fluidRow(
            column(12,
                   div(style = "background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
                        padding: 60px 40px;
                        border-radius: 10px;
                        margin: 50px 0;",
                       h2("How It Works",
                          style = "text-align: center; color: #667eea; margin-bottom: 50px; font-weight: 600; font-size: 2.2em;"),
                       fluidRow(
                         column(3,
                                div(style = "text-align: center;",
                                    div(style = "background: white;
                                         width: 90px;
                                         height: 90px;
                                         border-radius: 50%;
                                         margin: 0 auto 25px;
                                         display: flex;
                                         align-items: center;
                                         justify-content: center;
                                         box-shadow: 0 4px 10px rgba(0,0,0,0.15);",
                                        h2("1", style = "color: #667eea; margin: 0; font-size: 2.5em;")
                                    ),
                                    h4("Choose Method",
                                       style = "color: #333; margin-bottom: 15px; font-size: 1.3em;"),
                                    p("Select DPP, SAM, or Fisher's test",
                                      style = "color: #666; font-size: 1.05em;")
                                )
                         ),
                         column(3,
                                div(style = "text-align: center;",
                                    div(style = "background: white;
                                         width: 90px;
                                         height: 90px;
                                         border-radius: 50%;
                                         margin: 0 auto 25px;
                                         display: flex;
                                         align-items: center;
                                         justify-content: center;
                                         box-shadow: 0 4px 10px rgba(0,0,0,0.15);",
                                        h2("2", style = "color: #764ba2; margin: 0; font-size: 2.5em;")
                                    ),
                                    h4("Input Parameters",
                                       style = "color: #333; margin-bottom: 15px; font-size: 1.3em;"),
                                    p("Enter study design and historical data",
                                      style = "color: #666; font-size: 1.05em;")
                                )
                         ),
                         column(3,
                                div(style = "text-align: center;",
                                    div(style = "background: white;
                                         width: 90px;
                                         height: 90px;
                                         border-radius: 50%;
                                         margin: 0 auto 25px;
                                         display: flex;
                                         align-items: center;
                                         justify-content: center;
                                         box-shadow: 0 4px 10px rgba(0,0,0,0.15);",
                                        h2("3", style = "color: #f093fb; margin: 0; font-size: 2.5em;")
                                    ),
                                    h4("Run Analysis",
                                       style = "color: #333; margin-bottom: 15px; font-size: 1.3em;"),
                                    p("Calculate power and optimal design",
                                      style = "color: #666; font-size: 1.05em;")
                                )
                         ),
                         column(3,
                                div(style = "text-align: center;",
                                    div(style = "background: white;
                                         width: 90px;
                                         height: 90px;
                                         border-radius: 50%;
                                         margin: 0 auto 25px;
                                         display: flex;
                                         align-items: center;
                                         justify-content: center;
                                         box-shadow: 0 4px 10px rgba(0,0,0,0.15);",
                                        h2("4", style = "color: #4facfe; margin: 0; font-size: 2.5em;")
                                    ),
                                    h4("Interpret Results",
                                       style = "color: #333; margin-bottom: 15px; font-size: 1.3em;"),
                                    p("Get actionable insights for your trial",
                                      style = "color: #666; font-size: 1.05em;")
                                )
                         )
                       )
                   )
            )
          ),

          # Stats Section
          fluidRow(
            column(4,
                   div(style = "text-align: center; padding: 40px; background: white; border-radius: 10px; margin: 10px; box-shadow: 0 4px 15px rgba(0,0,0,0.1);",
                       h1("$1.2B", style = "color: #667eea; margin-bottom: 15px; font-size: 3.5em; font-weight: 700;"),
                       p("Average cost per drug", style = "color: #666; font-size: 1.2em;")
                   )
            ),
            column(4,
                   div(style = "text-align: center; padding: 40px; background: white; border-radius: 10px; margin: 10px; box-shadow: 0 4px 15px rgba(0,0,0,0.1);",
                       h1("10+ Years", style = "color: #764ba2; margin-bottom: 15px; font-size: 3.5em; font-weight: 700;"),
                       p("Development timeline", style = "color: #666; font-size: 1.2em;")
                   )
            ),
            column(4,
                   div(style = "text-align: center; padding: 40px; background: white; border-radius: 10px; margin: 10px; box-shadow: 0 4px 15px rgba(0,0,0,0.1);",
                       h1("Potential Time Savings", style = "color: #f093fb; margin-bottom: 15px; font-size: 3.5em; font-weight: 700;"),
                       p("Enhanced Statistical Inference", style = "color: #666; font-size: 1.2em;")
                   )
            )
          ),

          # Reference and Footer
          fluidRow(
            column(12,
                   div(style = "background: white;
                        padding: 35px;
                        border-radius: 10px;
                        box-shadow: 0 4px 15px rgba(0,0,0,0.1);
                        margin-top: 50px;",
                       h4(icon("book-open"), " Reference",
                          style = "color: #667eea; margin-bottom: 20px; font-size: 1.4em;"),
                       p(tags$a(
                         href = "https://onlinelibrary.wiley.com/doi/abs/10.1002/pst.2466",
                         target = "_blank",
                         "Lu Z, Toso J, Ayele G, He P. A Bayesian Hybrid Design With Borrowing From Historical Study. Pharm Stat. 2025 Mar-Apr;24(2):e2466. doi: 10.1002/pst.2466. Epub 2024 Dec 27. PMID: 39731333.",
                         style = "color: #667eea; text-decoration: none; font-size: 1.1em;"
                       )),
                       h4(icon("users"), " App Developer and Maintainer:",
                          style = "color: #667eea; margin-bottom: 20px; font-size: 1.4em;"),
                       p("Tanvi Mane, tm1049-AT-scarletmail.rutgers.edu. The App was created under IDSWG Oncology WG (https://oncologytrialdesign.org/). MIT License"),
                       hr(style = "border-color: #e0e0e0; margin: 25px 0;")
                   )
            )
          )
        ),

        # DPP Library Tab
        tabItem(
          tabName = "dpp_library",

          div(class = "content-header",
              h1(icon("book"), " Dynamic Power Prior - Parameter Library"),
              p("Click on any parameter below to view its definition and usage information.",
                style = "color: #666; font-size: 1.1em; margin-top: 10px;")
          ),

          fluidRow(
            column(4,
                   div(style = "background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       h4(icon("cog"), " General Parameters",
                          style = "color: #667eea; margin-bottom: 20px; font-size: 1.4em; font-weight: 600;"),
                       div(class = "term-link", actionLink("term_pt", HTML("<strong>pt</strong> - Experimental arm response rate"))),
                       div(class = "term-link", actionLink("term_nt", HTML("<strong>nt</strong> - Experimental arm sample size"))),
                       div(class = "term-link", actionLink("term_pc", HTML("<strong>pc</strong> - Control arm response rate"))),
                       div(class = "term-link", actionLink("term_nc", HTML("<strong>nc</strong> - Control arm sample size"))),
                       div(class = "term-link", actionLink("term_p_calib", HTML("<strong>pc.calib</strong> - Calibration response rate"))),
                       div(class = "term-link", actionLink("term_pch", HTML("<strong>pch</strong> - Historical control response rate"))),
                       div(class = "term-link", actionLink("term_nche", HTML("<strong>nche</strong> - Patients borrowed from historical"))),
                       div(class = "term-link", actionLink("term_nch", HTML("<strong>nch</strong> - Total historical control patients"))),
                       div(class = "term-link", actionLink("term_alpha", HTML("<strong>alpha</strong> - Type I error rate"))),
                       div(class = "term-link", actionLink("term_tau", HTML("<strong>tau</strong> - Significance threshold"))),
                       div(class = "term-link", actionLink("term_a0c", HTML("<strong>a0c</strong> - Control hyperprior alpha"))),
                       div(class = "term-link", actionLink("term_b0c", HTML("<strong>b0c</strong> - Control hyperprior beta"))),
                       div(class = "term-link", actionLink("term_a0t", HTML("<strong>a0t</strong> - Experimental hyperprior alpha"))),
                       div(class = "term-link", actionLink("term_b0t", HTML("<strong>b0t</strong> - Experimental hyperprior beta"))),
                       div(class = "term-link", actionLink("term_delta_threshold", HTML("<strong>delta_threshold</strong> - Borrowing threshold"))),
                       div(class = "term-link", actionLink("term_method", HTML("<strong>method</strong> - Dynamic borrowing method"))),
                       div(class = "term-link", actionLink("term_theta", HTML("<strong>theta</strong> - Method parameter (0-1)"))),
                       div(class = "term-link", actionLink("term_eta", HTML("<strong>eta</strong> - Method parameter (0-∞)"))),
                       div(class = "term-link", actionLink("term_datamat", HTML("<strong>datamat</strong> - Pre-simulated data matrix"))),
                       div(class = "term-link", actionLink("term_w0", HTML("<strong>w0</strong> - Prior power parameters"))),
                       div(class = "term-link", actionLink("term_nsim", HTML("<strong>nsim</strong> - Number of simulations"))),
                       div(class = "term-link", actionLink("term_seed", HTML("<strong>seed</strong> - Random seed")))
                   )
            ),
            column(4,
                   div(style = "background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       h4(icon("microscope"), " DPP Analysis Parameters",
                          style = "color: #11998e; margin-bottom: 20px; font-size: 1.4em; font-weight: 600;"),
                       div(class = "term-link", actionLink("dpp_analysis_term_w", HTML("<strong>w</strong> - Borrowing weight"))),
                       div(class = "term-link", actionLink("dpp_analysis_term_phat", HTML("<strong>phat_pt_larger_pc</strong> - Posterior probability"))),
                       div(class = "term-link", actionLink("dpp_analysis_term_apost_c_trial", HTML("<strong>apost_c_trial</strong> - Posterior alpha (trial)"))),
                       div(class = "term-link", actionLink("dpp_analysis_term_bpost_c_trial", HTML("<strong>bpost_c_trial</strong> - Posterior beta (trial)"))),
                       div(class = "term-link", actionLink("dpp_analysis_term_apost_c_hca", HTML("<strong>apost_c_hca</strong> - Posterior alpha (hybrid)"))),
                       div(class = "term-link", actionLink("dpp_analysis_term_bpost_c_hca", HTML("<strong>bpost_c_hca</strong> - Posterior beta (hybrid)"))),
                       div(class = "term-link", actionLink("dpp_analysis_term_apost_t", HTML("<strong>apost_t</strong> - Posterior alpha (experimental)"))),
                       div(class = "term-link", actionLink("dpp_analysis_term_bpost_t", HTML("<strong>bpost_t</strong> - Posterior beta (experimental)")))
                   )
            ),
            column(4,
                   div(class = "definition-box",
                       h4(icon("info-circle"), " Parameter Definition",
                          style = "color: #667eea; margin-bottom: 20px; font-size: 1.4em; font-weight: 600;"),
                       uiOutput("definition_output"),
                       div(style = "margin-top: 20px; padding: 15px; background: #f8f9fc; border-radius: 6px; border-left: 4px solid #f093fb;",
                           p(icon("lightbulb"), strong(" Tip: "), "Click on any parameter name to see its detailed definition here.",
                             style = "margin: 0; color: #666;")
                       )
                   )
            )
          )
        ),

        # DPP Design Tab
        tabItem(
          tabName = "dpp_design",

          div(class = "content-header",
              h1(icon("calculator"), " Single Study Design - Dynamic Power Prior"),
              p("Configure your study parameters below to calculate power, tau, and other design metrics.",
                style = "color: #666; font-size: 1.1em; margin-top: 10px;")
          ),

          fluidRow(
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("flask"), " Current Study Parameters",
                          style = "color: #667eea; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       fluidRow(
                         column(6, numericInput("nt", "Experimental Arm Sample Size (nt)", value = 50, min = 1)),
                         column(6, numericInput("pt", "Experimental Arm Response Rate (pt)", value = 0.5, min = 0, max = 1, step = 0.01))
                       ),
                       fluidRow(
                         column(6, numericInput("nc", "Control Arm Sample Size (nc)", value = 50, min = 1)),
                         column(6, numericInput("pc", "Control Arm Response Rate (pc)", value = 0.3, min = 0, max = 1, step = 0.01))
                       )
                   ),

                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("database"), " Historical Control Data",
                          style = "color: #764ba2; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       fluidRow(
                         column(6, numericInput("nch", "Historical Control Sample Size (nch)", value = 100, min = 1)),
                         column(6, numericInput("pch", "Historical Control Response Rate (pch)", value = 0.3, min = 0, max = 1, step = 0.01))
                       )
                   ),

                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("exchange-alt"), " Bayesian Borrowing",
                          style = "color: #11998e; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("nche", "Max Patients from Historical Study (nche)", value = 50, min = 1),
                       numericInput("delta_threshold", "Delta Threshold", value = 0.1, min = 0, max = 1, step = 0.01),
                       selectInput("method", "Dynamic Control Method",
                                   choices = c("Empirical Bayes", "Bayesian p", "Generalized BC", "JSD")),
                       uiOutput("theta_ui"),
                       uiOutput("eta_ui")
                   )
            ),

            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("chart-bar"), " Statistical Parameters",
                          style = "color: #667eea; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("alpha", "One-Sided Type I Error", value = 0.1, min = 0.001, max = 0.5),
                       numericInput("p_calib", "Response Rate for Calibration (pc.calib)", value = 0.3, min = 0, max = 1, step = 0.01)
                   ),

                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("sliders-h"), " Hyperpriors (Beta Distribution)",
                          style = "color: #764ba2; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       h5("Control Arm", style = "color: #666; margin-top: 15px;"),
                       fluidRow(
                         column(6, numericInput("a0c", "a0c", value = 0.001, min = 0)),
                         column(6, numericInput("b0c", "b0c", value = 0.001, min = 0))
                       ),
                       h5("Experimental Arm", style = "color: #666; margin-top: 15px;"),
                       fluidRow(
                         column(6, numericInput("a0t", "a0t", value = 0.001, min = 0)),
                         column(6, numericInput("b0t", "b0t", value = 0.001, min = 0))
                       )
                   ),

                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("random"), " Simulation Settings",
                          style = "color: #11998e; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       fluidRow(
                         column(6, numericInput("nsim", "Number of Simulations", value = 100000, min = 100, step = 1000)),
                         column(6, textInput("seed", "Random Seed", value = "2000"))
                       )
                   )
            )
          ),

          fluidRow(
            column(12,
                   div(style = "text-align: center; margin: 30px 0;",
                       actionButton("run_dpp", "Run Study Design Analysis",
                                    class = "btn-primary btn-lg",
                                    icon = icon("play-circle"),
                                    style = "padding: 15px 50px; font-size: 1.3em;")
                   ),
                   div(id = "analysis_status_message", style = "text-align: center; margin-top: 20px;")
            )
          ),

          # ===========================================================================
          # RESULTS SECTION - UPDATED WITH BOTH METRICS
          # ===========================================================================

          # Key Metrics Value Boxes
          fluidRow(
            column(4, valueBoxOutput("dppPowerBox", width = NULL)),
            column(4, valueBoxOutput("dppTauBox", width = NULL)),
            column(4, valueBoxOutput("dppDeltaBoundBox", width = NULL))
          ),

          # PMD Statistics Section
          fluidRow(
            column(12,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin: 20px 0;",
                       h4(icon("chart-line"), " Posterior Mean Difference (PMD) Statistics",
                          style = "color: #667eea; font-weight: 600; margin-bottom: 25px;"),

                       fluidRow(
                         column(4,
                                div(style = "text-align: center; padding: 20px; border-right: 1px solid #e8e9ef;",
                                    h5("PMD Mean", style = "color: #764ba2; font-weight: 600; margin-bottom: 10px;"),
                                    div(style = "font-size: 2.5em; color: #764ba2; font-weight: 700;",
                                        textOutput("dpp_pmd_mean")),
                                    p("Mean difference between hybrid control and concurrent control",
                                      style = "color: #666; margin-top: 10px; font-size: 0.9em;")
                                )
                         ),
                         column(4,
                                div(style = "text-align: center; padding: 20px; border-right: 1px solid #e8e9ef;",
                                    h5("PMD Standard Deviation", style = "color: #f093fb; font-weight: 600; margin-bottom: 10px;"),
                                    div(style = "font-size: 2.5em; color: #f093fb; font-weight: 700;",
                                        textOutput("dpp_pmd_sd")),
                                    p("Variability of the posterior mean difference",
                                      style = "color: #666; margin-top: 10px; font-size: 0.9em;")
                                )
                         ),
                         column(4,
                                div(style = "text-align: center; padding: 20px;",
                                    h5("95% Credible Interval", style = "color: #11998e; font-weight: 600; margin-bottom: 10px;"),
                                    div(style = "font-size: 1.8em; color: #11998e; font-weight: 700;",
                                        textOutput("dpp_pmd_ci")),
                                    p("95% credible interval for PMD",
                                      style = "color: #666; margin-top: 10px; font-size: 0.9em;")
                                )
                         )
                       )
                   )
            )
          ),

          # PMD Distribution Plot
          fluidRow(
            column(12,
                   div(style = "background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       h4(icon("chart-bar"), " Posterior Mean Difference (PMD) Distribution",
                          style = "color: #f093fb; font-weight: 600; margin-bottom: 20px;"),
                       p("This plot shows the density of the posterior mean difference between the hybrid control arm and the concurrent control arm across simulated trials.",
                         style = "color: #666; margin-bottom: 20px;"),
                       plotOutput("dpp_plot_pmd", height = "500px")
                   )
            )
          )
        ),

        # DPP Comparative Analysis Tab
        tabItem(
          tabName = "dpp_table",

          div(class = "content-header",
              h1(icon("table"), " Comparative Analysis - Multiple Scenarios"),
              p("Compare Type I error, Power, and PMD across different control response rates and historical borrowing amounts.",
                style = "color: #666; font-size: 1.1em; margin-top: 10px;")
          ),

          fluidRow(
            column(12,
                   div(style = "background: #fff3cd; padding: 20px; border-radius: 8px; border-left: 5px solid #ffc107; margin-bottom: 20px;",
                       h5(icon("exclamation-triangle"), strong(" Important Note"), style = "color: #856404; margin-bottom: 10px;"),
                       p(strong("Note:"), " This analysis may take several minutes depending on the number of scenarios and simulations.",
                         style = "color: #856404; margin: 0; font-size: 1.05em;")
                   )
            )
          ),

          fluidRow(
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("flask"), " Study Design Parameters",
                          style = "color: #667eea; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("table_nt", "Experimental Arm Sample Size", value = 45, min = 1),
                       numericInput("table_nc", "Control Arm Sample Size", value = 45, min = 1),
                       textInput("table_pc_values", "Control Arm Response Rates (comma-separated)",
                                 value = "0.15,0.2,0.25,0.3,0.35,0.4,0.45",
                                 placeholder = "e.g., 0.15,0.2,0.25"),
                       numericInput("table_pch", "Historical Control Response Rate", value = 0.3, min = 0, max = 1, step = 0.01),
                       numericInput("table_nch", "Historical Control Sample Size", value = 180, min = 1)
                   )
            ),
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("exchange-alt"), " Borrowing & Analysis Parameters",
                          style = "color: #764ba2; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       textInput("table_nche_values", "nche Values (comma-separated)",
                                 value = "45,90,135,180",
                                 placeholder = "e.g., 45,90,135,180"),
                       numericInput("table_delta_threshold", "Delta Threshold", value = 0.1, min = 0, max = 1, step = 0.01),
                       numericInput("table_alpha", "Type I Error", value = 0.1, min = 0.001, max = 0.5)
                   ),
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("sliders-h"), " Advanced Settings",
                          style = "color: #11998e; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       fluidRow(
                         column(6, numericInput("table_a0c", "Control a0c", value = 0.001, min = 0)),
                         column(6, numericInput("table_b0c", "Control b0c", value = 0.001, min = 0))
                       ),
                       fluidRow(
                         column(6, numericInput("table_a0t", "Experimental a0t", value = 0.001, min = 0)),
                         column(6, numericInput("table_b0t", "Experimental b0t", value = 0.001, min = 0))
                       ),
                       numericInput("table_nsim", "Number of Simulations", value = 10000, min = 100, step = 1000)
                   )
            )
          ),

          fluidRow(
            column(12,
                   div(style = "text-align: center; margin: 30px 0;",
                       actionButton("run_dpp_table", "Run Comparative Analysis",
                                    class = "btn-warning btn-lg",
                                    icon = icon("table"),
                                    style = "padding: 15px 50px; font-size: 1.3em;")
                   ),
                   div(id = "dpp_table_status_message", style = "text-align: center; margin-top: 20px;")
            )
          ),

          fluidRow(
            column(12,
                   div(style = "background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       tabBox(
                         title = "Results",
                         width = NULL,
                         tabPanel("Type I Error",
                                  p("Rows: Control response rates (pc) | Columns: Borrowing amounts (nche)"),
                                  DT::dataTableOutput("dpp_table_typeI")
                         ),
                         tabPanel("Power",
                                  p("Power assuming treatment effect = control rate + 0.2"),
                                  DT::dataTableOutput("dpp_table_power")
                         ),
                         tabPanel("Posterior Mean Difference",
                                  p("Mean difference between hybrid and concurrent control"),
                                  DT::dataTableOutput("dpp_table_pmd")
                         ),
                         tabPanel("SD of PMD",
                                  DT::dataTableOutput("dpp_table_sd_pmd")
                         )
                       )
                   )
            )
          )
        ),

        # DPP Statistical Analysis Tab
        tabItem(
          tabName = "dpp_analysis",

          div(class = "content-header",
              h1(icon("microscope"), " Statistical Analysis - Analyze Observed Trial Data"),
              p("Enter your observed trial data to perform DPP analysis and get posterior distributions.",
                style = "color: #666; font-size: 1.1em; margin-top: 10px;")
          ),

          fluidRow(
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("user-friends"), " Observed Trial Data",
                          style = "color: #667eea; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       fluidRow(
                         column(6, numericInput("dpp_analysis_rt", "Responders in Experimental Arm (rt)", value = 30, min = 0)),
                         column(6, numericInput("dpp_analysis_nt", "Experimental Arm Sample Size (nt)", value = 41, min = 1))
                       ),
                       fluidRow(
                         column(6, numericInput("dpp_analysis_rc", "Responders in Control Arm (rc)", value = 15, min = 0)),
                         column(6, numericInput("dpp_analysis_nc", "Control Arm Sample Size (nc)", value = 44, min = 1))
                       ),
                       div(style = "margin-top: 15px; padding: 15px; background-color: #f8f9fc; border-radius: 6px; border-left: 4px solid #667eea;",
                           h5("Calculated Response Rates", style = "color: #667eea; margin-bottom: 10px;"),
                           fluidRow(
                             column(6,
                                    strong("Experimental (pt):"),
                                    verbatimTextOutput("dpp_analysis_pt_display", placeholder = TRUE)
                             ),
                             column(6,
                                    strong("Control (pc):"),
                                    verbatimTextOutput("dpp_analysis_pc_display", placeholder = TRUE)
                             )
                           )
                       )
                   ),

                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("database"), " Historical Control Data",
                          style = "color: #764ba2; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       fluidRow(
                         column(6, numericInput("dpp_analysis_nch", "Historical Control Sample Size (nch)", value = 87, min = 1)),
                         column(6, numericInput("dpp_analysis_pch", "Historical Control Response Rate (pch)", value = 0.3, min = 0, max = 1, step = 0.01))
                       )
                   )
            ),

            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("exchange-alt"), " Bayesian Borrowing Configuration",
                          style = "color: #11998e; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("dpp_analysis_nche", "Equivalent Patients Borrowed (nche)", value = 41, min = 1),
                       numericInput("dpp_analysis_delta_threshold", "Delta Threshold", value = 0.1, min = 0, max = 0.5, step = 0.01),
                       selectInput("dpp_analysis_method", "Dynamic Control Method",
                                   choices = c("Empirical Bayes", "Bayesian p", "Generalized BC", "JSD")),
                       uiOutput("dpp_analysis_theta_ui"),
                       uiOutput("dpp_analysis_eta_ui")
                   ),

                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("sliders-h"), " Hyperpriors",
                          style = "color: #667eea; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       fluidRow(
                         column(6, numericInput("dpp_analysis_a0c", "Control a0c", value = 0.001, min = 0)),
                         column(6, numericInput("dpp_analysis_b0c", "Control b0c", value = 0.001, min = 0))
                       ),
                       fluidRow(
                         column(6, numericInput("dpp_analysis_a0t", "Experimental a0t", value = 0.001, min = 0)),
                         column(6, numericInput("dpp_analysis_b0t", "Experimental b0t", value = 0.001, min = 0))
                       )
                   )
            )
          ),

          fluidRow(
            column(12,
                   div(style = "text-align: center; margin: 30px 0;",
                       actionButton("run_dpp_analysis", "Run DPP Analysis",
                                    class = "btn-success btn-lg",
                                    icon = icon("microscope"),
                                    style = "padding: 15px 50px; font-size: 1.3em;")
                   ),
                   div(id = "dpp_analysis_status_message", style = "text-align: center; margin-top: 20px;")
            )
          ),

          fluidRow(
            column(12,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("chart-bar"), " Analysis Results",
                          style = "color: #667eea; font-weight: 600; margin-bottom: 20px;"),
                       verbatimTextOutput("dppAnalysisResult"),
                       div(style = "margin-top: 15px; padding: 15px; background: #f8f9fc; border-radius: 6px; border-left: 4px solid #f093fb;",
                           p(icon("lightbulb"), strong(" Note:"), " For parameter definitions, check the Library tab.",
                             style = "margin: 0; color: #666;")
                       )
                   )
            )
          ),

          fluidRow(
            column(12,
                   div(style = "background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("chart-area"), " Posterior Distributions",
                          style = "color: #764ba2; font-weight: 600; margin-bottom: 20px;"),
                       plotOutput("plotDPP", height = "500px")
                   )
            )
          ),

          fluidRow(
            column(4,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("vial"), " Concurrent Control",
                          style = "color: #667eea; font-weight: 600; margin-bottom: 15px;"),
                       verbatimTextOutput("concurrent_summary")
                   )
            ),
            column(4,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("layer-group"), " Hybrid Control",
                          style = "color: #764ba2; font-weight: 600; margin-bottom: 15px;"),
                       verbatimTextOutput("hybrid_control_summary")
                   )
            ),
            column(4,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("flask"), " Experimental Treatment",
                          style = "color: #11998e; font-weight: 600; margin-bottom: 15px;"),
                       verbatimTextOutput("experimental_summary")
                   )
            )
          ),

          fluidRow(
            column(12,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       h4(icon("check-circle"), " Statistical Conclusion",
                          style = "color: #667eea; font-weight: 600; margin-bottom: 15px;"),
                       verbatimTextOutput("statistical_conclusion")
                   )
            )
          )
        ),

        # SAM Library Tab
        tabItem(
          tabName = "sam_library",

          div(class = "content-header",
              h1(icon("book"), " SAM Prior - Parameter Library"),
              p("Click on any parameter below to view its definition and usage information.",
                style = "color: #666; font-size: 1.1em; margin-top: 10px;")
          ),

          fluidRow(
            column(4,
                   div(style = "background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       h4(icon("chart-bar"), " SAM Prior Plot Parameters",
                          style = "color: #667eea; margin-bottom: 20px; font-size: 1.4em; font-weight: 600;"),
                       div(class = "term-link", actionLink("sam_plot_library_term_alpha_hist", HTML("<strong>alpha_hist</strong> - Historical prior alpha"))),
                       div(class = "term-link", actionLink("sam_plot_library_term_beta_hist", HTML("<strong>beta_hist</strong> - Historical prior beta"))),
                       div(class = "term-link", actionLink("sam_plot_library_term_n_control", HTML("<strong>n_control</strong> - Control arm sample size"))),
                       div(class = "term-link", actionLink("sam_plot_library_term_p_control", HTML("<strong>p_control</strong> - Simulated control rate"))),
                       div(class = "term-link", actionLink("sam_plot_library_term_nf_alpha", HTML("<strong>nf_alpha</strong> - Non-informative prior alpha"))),
                       div(class = "term-link", actionLink("sam_plot_library_term_nf_beta", HTML("<strong>nf_beta</strong> - Non-informative prior beta"))),
                       div(class = "term-link", actionLink("sam_plot_library_term_delta", HTML("<strong>delta</strong> - Clinically meaningful difference")))
                   )
            ),
            column(4,
                   div(style = "background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       h4(icon("balance-scale"), " SAM Weight Parameters",
                          style = "color: #11998e; margin-bottom: 20px; font-size: 1.4em; font-weight: 600;"),
                       div(class = "term-link", actionLink("sam_weight_library_term_alpha_hist_w", HTML("<strong>alpha_hist</strong> - Historical prior alpha"))),
                       div(class = "term-link", actionLink("sam_weight_library_term_beta_hist_w", HTML("<strong>beta_hist</strong> - Historical prior beta"))),
                       div(class = "term-link", actionLink("sam_weight_library_term_n_control_w", HTML("<strong>n_control</strong> - Control arm sample size"))),
                       div(class = "term-link", actionLink("sam_weight_library_term_r_control_w", HTML("<strong>r_control</strong> - Number of responses"))),
                       div(class = "term-link", actionLink("sam_weight_library_term_delta_w", HTML("<strong>delta</strong> - Clinically meaningful difference"))),
                       div(class = "term-link", actionLink("sam_weight_library_term_method_w", HTML("<strong>method</strong> - Weight calculation method"))),
                       div(class = "term-link", actionLink("sam_weight_library_term_prior_odds", HTML("<strong>prior_odds</strong> - Prior odds (PPR only)")))
                   )
            ),
            column(4,
                   div(class = "definition-box",
                       h4(icon("info-circle"), " Parameter Definition",
                          style = "color: #667eea; margin-bottom: 20px; font-size: 1.4em; font-weight: 600;"),
                       uiOutput("sam_plot_definition_output"),
                       div(style = "margin-top: 20px; padding: 15px; background: #f8f9fc; border-radius: 6px; border-left: 4px solid #f093fb;",
                           p(icon("lightbulb"), strong(" Tip: "), "Click on any parameter name to see its detailed definition here.",
                             style = "margin: 0; color: #666;")
                       )
                   )
            )
          )
        ),

        # SAM Design Tab
        tabItem(
          tabName = "sam_design",

          div(class = "content-header",
              h1(icon("calculator"), " Single Study Design - SAM Prior"),
              p("Configure your study parameters below to generate SAM prior plots and calculate weights.",
                style = "color: #666; font-size: 1.1em; margin-top: 10px;")
          ),

          fluidRow(
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("user-friends"), " Current Study Parameters",
                          style = "color: #667eea; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       fluidRow(
                         column(6, numericInput("sam_nt", "Experimental Arm Sample Size (nt)", value = 50, min = 1, step = 1)),
                         column(6, numericInput("sam_pt", "Experimental Arm Response Rate (pt)", value = 0.5, min = 0, max = 1, step = 0.01))
                       ),
                       fluidRow(
                         column(6, numericInput("sam_nc", "Control Arm Sample Size (nc)", value = 50, min = 1, step = 1)),
                         column(6, numericInput("sam_pc", "Control Arm Response Rate (pc)", value = 0.3, min = 0, max = 1, step = 0.01))
                       )
                   ),

                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("database"), " Historical Control Data",
                          style = "color: #764ba2; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       fluidRow(
                         column(6, numericInput("sam_nch", "Historical Control Sample Size (nch)", value = 100, min = 0, step = 1)),
                         column(6, numericInput("sam_pch", "Historical Control Response Rate (pch)", value = 0.3, min = 0, max = 1, step = 0.01))
                       )
                   ),

                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("chart-bar"), " SAM Prior Plot Settings",
                          style = "color: #11998e; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("sam_sim_control_rate", "Simulated Control Rate", value = 0.4, min = 0, max = 1, step = 0.01),
                       h5("Non-Informative Prior", style = "color: #666; margin-top: 15px;"),
                       fluidRow(
                         column(6, numericInput("sam_alpha_noninf", "Alpha (a)", value = 1, min = 0, step = 1)),
                         column(6, numericInput("sam_beta_noninf", "Beta (b)", value = 1, min = 0, step = 1))
                       )
                   )
            ),

            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("sliders-h"), " Design Parameters",
                          style = "color: #667eea; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("sam_alpha", "One-Sided Type I Error", value = 0.1, min = 0.001, max = 0.5, step = 0.001),
                       numericInput("sam_pc_calib", "Response Rate for Calibration (pc.calib)", value = 0.3, min = 0, max = 1, step = 0.01),
                       numericInput("sam_delta", HTML("Clinically Meaningful Difference (&delta;)"), value = 0.15, min = 0, max = 1, step = 0.01)
                   ),

                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("random"), " Simulation Settings",
                          style = "color: #764ba2; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       fluidRow(
                         column(6, numericInput("sam_n_sim", "Number of Simulations", value = 100000, min = 1000, step = 1000)),
                         column(6, numericInput("sam_seed", "Random Seed", value = 2000, min = 1, step = 1))
                       )
                   )
            )
          ),

          fluidRow(
            column(12,
                   div(style = "text-align: center; margin: 30px 0;",
                       actionButton("generate_sam_prior_plot", "Generate SAM Prior Plot & Weight",
                                    class = "btn-success btn-lg",
                                    icon = icon("chart-line"),
                                    style = "padding: 15px 50px; font-size: 1.3em;")
                   )
            )
          ),

          fluidRow(
            column(6,
                   div(style = "background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
                      padding: 30px;
                      border-radius: 12px;
                      box-shadow: 0 8px 20px rgba(102, 126, 234, 0.3);
                      color: white;
                      min-height: 220px;",
                       div(style = "display: flex; align-items: center; margin-bottom: 20px;",
                           div(style = "background: rgba(255,255,255,0.2);
                              padding: 15px;
                              border-radius: 50%;
                              margin-right: 15px;
                              width: 60px;
                              height: 60px;
                              display: flex;
                              align-items: center;
                              justify-content: center;",
                               icon("calculator", style = "font-size: 28px;")
                           ),
                           h4("Calculated Beta Prior Parameters",
                              style = "margin: 0; font-weight: 600; font-size: 1.4em;")
                       ),
                       p("Informative prior derived from historical control data",
                         style = "opacity: 0.9; margin-bottom: 20px; font-size: 1.05em;"),
                       div(style = "background: rgba(255,255,255,0.15);
                          padding: 20px;
                          border-radius: 8px;
                          border-left: 4px solid rgba(255,255,255,0.5);",
                           uiOutput("sam_calculated_prior")
                       )
                   )
            ),

            column(6,
                   div(style = "background: linear-gradient(135deg, #11998e 0%, #38ef7d 100%);
                      padding: 30px;
                      border-radius: 12px;
                      box-shadow: 0 8px 20px rgba(17, 153, 142, 0.3);
                      color: white;
                      min-height: 220px;",
                       div(style = "display: flex; align-items: center; margin-bottom: 20px;",
                           div(style = "background: rgba(255,255,255,0.2);
                              padding: 15px;
                              border-radius: 50%;
                              margin-right: 15px;
                              width: 60px;
                              height: 60px;
                              display: flex;
                              align-items: center;
                              justify-content: center;",
                               icon("balance-scale", style = "font-size: 28px;")
                           ),
                           h4("SAM Weight",
                              style = "margin: 0; font-weight: 600; font-size: 1.4em;")
                       ),
                       p("Optimal mixture of informative and non-informative priors",
                         style = "opacity: 0.9; margin-bottom: 20px; font-size: 1.05em;"),
                       div(style = "background: rgba(255,255,255,0.15);
                          padding: 20px;
                          border-radius: 8px;
                          border-left: 4px solid rgba(255,255,255,0.5);",
                           uiOutput("sam_weight_display")
                       )
                   )
            )
          ),

          fluidRow(
            column(12,
                   div(style = "background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-top: 20px;",
                       h4(icon("chart-area"), " SAM Prior Distribution",
                          style = "color: #f093fb; font-weight: 600; margin-bottom: 20px;"),
                       plotOutput("sam_prior_plot", height = "500px")
                   )
            )
          )
        ),

        # SAM Comparative Analysis Tab
        tabItem(
          tabName = "sam_table",

          div(class = "content-header",
              h1(icon("table"), " SAM Prior Comparative Analysis"),
              p("Compare SAMprior, MAP, and Non-informative methods across different control response rates using calibrated thresholds.",
                style = "color: #666; font-size: 1.1em; margin-top: 10px;")
          ),

          fluidRow(
            column(12,
                   div(style = "background: #fff3cd; padding: 20px; border-radius: 8px; border-left: 5px solid #ffc107; margin-bottom: 20px;",
                       h5(icon("exclamation-triangle"), strong(" Important Note"), style = "color: #856404; margin-bottom: 10px;"),
                       p("SAM analysis is computationally intensive. The new calibrated approach performs two-step simulations:",
                         style = "color: #856404; margin: 0; font-size: 1.05em;"),
                       tags$ul(
                         tags$li("Calibration step: Determines thresholds under the null hypothesis"),
                         tags$li("Main simulation: Evaluates operating characteristics using calibrated thresholds"),
                         style = "color: #856404; margin-top: 10px;"
                       ),
                       p(strong("Recommended:"), " Start with 5,000-10,000 simulations for testing, then increase to 100,000 for final analysis.",
                         style = "color: #856404; margin-top: 10px; font-size: 1.05em;")
                   )
            )
          ),

          fluidRow(
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       h4(icon("flask"), " Study Design Parameters",
                          style = "color: #667eea; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("sam_table_nt", "Experimental Arm Sample Size (nt)", value = 45, min = 1),
                       numericInput("sam_table_nc", "Control Arm Sample Size (nc)", value = 45, min = 1),
                       textInput("sam_table_pc_values", "Control Response Rates (comma-separated)",
                                 value = "0.15,0.2,0.25,0.3,0.35,0.4,0.45",
                                 placeholder = "e.g., 0.15,0.2,0.25"),
                       helpText("These values will be used for both calibration and main simulations"),
                       numericInput("sam_table_pch", "Historical Control Response Rate (pch)", value = 0.3, min = 0, max = 1, step = 0.01),
                       numericInput("sam_table_nche", "Historical Control Sample Size (nche)", value = 180, min = 1),
                       helpText("This represents the equivalent number of patients borrowed from historical data")
                   )
            ),
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       h4(icon("sliders-h"), " SAM Parameters",
                          style = "color: #764ba2; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("sam_table_delta", "Delta Threshold (CSD)", value = 0.1, min = 0, max = 1, step = 0.01),
                       helpText("Clinically Significant Difference threshold for SAM prior"),
                       numericInput("sam_table_typeIER", "Target Type I Error Rate", value = 0.1, min = 0.001, max = 0.5, step = 0.01),
                       helpText("Target Type I Error to control during calibration (default: 0.10)"),
                       numericInput("sam_table_nsim", "Number of Simulations", value = 10000, min = 100, step = 1000),
                       helpText("Used for both calibration and main simulation steps"),
                       div(style = "margin-top: 15px; padding: 15px; background: #e7f3ff; border-radius: 6px; border-left: 4px solid #2196F3;",
                           p(icon("info-circle"), strong(" Calibration Process: "),
                             "The function will first calibrate thresholds under the null hypothesis, then use these to evaluate Type I Error and Power.",
                             style = "margin: 0; color: #0d47a1; font-size: 0.95em;")
                       )
                   )
            )
          ),

          fluidRow(
            column(12,
                   div(style = "text-align: center; margin: 30px 0;",
                       actionButton("run_sam_table", "Run SAM Comparative Analysis (Calibrated)",
                                    class = "btn-warning btn-lg",
                                    icon = icon("table"),
                                    style = "padding: 15px 50px; font-size: 1.3em;")
                   ),
                   div(id = "sam_table_status_message", style = "text-align: center; margin-top: 20px;")
            )
          ),

          fluidRow(
            column(12,
                   div(style = "background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       tabBox(
                         title = "SAM Analysis Results (Calibrated)",
                         width = NULL,
                         tabPanel(
                           title = tagList(icon("exclamation-triangle"), " Type I Error"),
                           p("Type I error rates for SAMprior, MAP, and Non-informative methods using calibrated thresholds",
                             style = "color: #666; margin-bottom: 15px;"),
                           helpText("Note: These should be close to the target Type I Error Rate specified above"),
                           DT::dataTableOutput("sam_table_typeI")
                         ),
                         tabPanel(
                           title = tagList(icon("bolt"), " Power"),
                           p("Power assuming treatment effect = control rate + 0.2, using calibrated thresholds",
                             style = "color: #666; margin-bottom: 15px;"),
                           helpText("Rows: Control response rates | Columns: SAMprior, MAP, Non-informative methods"),
                           DT::dataTableOutput("sam_table_power")
                         ),
                         tabPanel(
                           title = tagList(icon("chart-line"), " Posterior Mean Difference"),
                           p("Mean difference between methods and non-informative approach",
                             style = "color: #666; margin-bottom: 15px;"),
                           helpText("This shows the bias introduced by borrowing from historical data"),
                           DT::dataTableOutput("sam_table_pmd")
                         ),
                         tabPanel(
                           title = tagList(icon("signal"), " SD of PMD"),
                           p("Standard deviation of posterior mean difference",
                             style = "color: #666; margin-bottom: 15px;"),
                           helpText("This measures the variability of the bias across simulations"),
                           DT::dataTableOutput("sam_table_sd_pmd")
                         )
                       )
                   )
            )
          )
        ),



        # Fisher Library Tab
        tabItem(
          tabName = "fisher_library",

          div(class = "content-header",
              h1(icon("book"), " Fisher's Exact Method - Parameter Library"),
              p("Click on any parameter below to view its definition and usage information.",
                style = "color: #666; font-size: 1.1em; margin-top: 10px;")
          ),

          fluidRow(
            column(6,
                   div(style = "background: white; padding: 30px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       h4(icon("dice"), " Fisher's Exact Test Parameters",
                          style = "color: #667eea; margin-bottom: 20px; font-size: 1.4em; font-weight: 600;"),
                       div(class = "term-link", actionLink("fisher_term_Yc", "Yc (fisher)")),
                       div(class = "term-link", actionLink("fisher_term_nc_fisher", "nc (fisher)")),
                       div(class = "term-link", actionLink("fisher_term_Yt", "Yt (fisher)")),
                       div(class = "term-link", actionLink("fisher_term_nt_fisher", "nt (fisher)")),
                       div(class = "term-link", actionLink("fisher_term_alternative", "alternative (fisher)")),
                       div(class = "term-link", actionLink("fisher_term_pc_bound", "pc (fisher.bound)")),
                       div(class = "term-link", actionLink("fisher_term_nc_bound", "nc (fisher.bound)")),
                       div(class = "term-link", actionLink("fisher_term_nt_bound", "nt (fisher.bound)")),
                       div(class = "term-link", actionLink("fisher_term_alpha_bound", "alpha (fisher.bound)")),
                       div(class = "term-link", actionLink("fisher_term_pt_power", "pt (fisher.power)")),
                       div(class = "term-link", actionLink("fisher_term_nt_power", "nt (fisher.power)")),
                       div(class = "term-link", actionLink("fisher_term_pc_power", "pc (fisher.power)")),
                       div(class = "term-link", actionLink("fisher_term_nc_power", "nc (fisher.power)")),
                       div(class = "term-link", actionLink("fisher_term_alpha_power", "alpha (fisher.power)")),
                       div(class = "term-link", actionLink("fisher_term_nsim_power", "nsim (fisher.power)")),
                       div(class = "term-link", actionLink("fisher_term_seed_power", "seed (fisher.power)"))
                   )
            ),
            column(6,
                   div(class = "definition-box",
                       h4(icon("info-circle"), " Parameter Definition",
                          style = "color: #667eea; margin-bottom: 20px; font-size: 1.4em; font-weight: 600;"),
                       uiOutput("fisher_definition_output"),
                       div(style = "margin-top: 20px; padding: 15px; background: #f8f9fc; border-radius: 6px; border-left: 4px solid #f093fb;",
                           p(icon("lightbulb"), strong(" Tip: "), "Click on any parameter name to see its detailed definition here.",
                             style = "margin: 0; color: #666;")
                       )
                   )
            )
          )
        ),

        # Fisher Power Tab
        tabItem(
          tabName = "fisher_power",

          div(class = "content-header",
              h1(icon("bolt"), " Fisher's Exact Test - Power Analysis"),
              p("Calculate statistical power using Fisher's Exact test for your study design.",
                style = "color: #666; font-size: 1.1em; margin-top: 10px;")
          ),

          fluidRow(
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("flask"), " Study Parameters",
                          style = "color: #667eea; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("fp_pt", "Experimental Arm Response Rate (pt)", value = 0.5, min = 0, max = 1, step = 0.01),
                       numericInput("fp_nt", "Experimental Arm Sample Size (nt)", value = 40, min = 1),
                       numericInput("fp_pc", "Control Arm Response Rate (pc)", value = 0.3, min = 0, max = 1, step = 0.01),
                       numericInput("fp_nc", "Control Arm Sample Size (nc)", value = 40, min = 1)
                   )
            ),
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("sliders-h"), " Analysis Parameters",
                          style = "color: #764ba2; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("fp_alpha", "One-Sided Type I Error (α)", value = 0.1, min = 0.001, max = 0.5),
                       numericInput("fp_nsim", "Number of Simulations", value = 100000, min = 100, step = 1000),
                       numericInput("fp_seed", "Random Seed", value = 2000, min = 1)
                   )
            )
          ),

          fluidRow(
            column(12,
                   div(style = "text-align: center; margin: 30px 0;",
                       actionButton("run_fisher_power", "Calculate Fisher Power",
                                    class = "btn-primary btn-lg",
                                    icon = icon("bolt"),
                                    style = "padding: 15px 50px; font-size: 1.3em;")
                   ),
                   div(id = "fisher_power_status_message", style = "text-align: center; margin-top: 20px;")
            )
          ),

          fluidRow(
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("chart-line"), " Power Result",
                          style = "color: #11998e; font-weight: 600; margin-bottom: 15px;"),
                       verbatimTextOutput("fisherPowerResult")
                   )
            ),
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("lightbulb"), " Interpretation",
                          style = "color: #667eea; font-weight: 600; margin-bottom: 15px;"),
                       htmlOutput("fisherPowerConclusion")
                   )
            )
          )
        ),

        # Fisher Bound Tab
        tabItem(
          tabName = "fisher_bound",

          div(class = "content-header",
              h1(icon("border-all"), " Fisher's Exact Test - Boundary Analysis"),
              p("Determine the minimum number of responders needed for statistical significance.",
                style = "color: #666; font-size: 1.1em; margin-top: 10px;")
          ),

          fluidRow(
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("flask"), " Study Parameters",
                          style = "color: #667eea; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("fb_pc", "Control Arm Response Rate (pc)", value = 0.3, min = 0, max = 1, step = 0.01),
                       numericInput("fb_nc", "Control Arm Sample Size (nc)", value = 40, min = 1)
                   )
            ),
            column(6,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("sliders-h"), " Analysis Parameters",
                          style = "color: #764ba2; margin-bottom: 20px; font-weight: 600; border-bottom: 2px solid #f8f9fc; padding-bottom: 10px;"),
                       numericInput("fb_nt", "Experimental Arm Sample Size (nt)", value = 40, min = 1),
                       numericInput("fb_alpha", "Significance Threshold (α)", value = 0.1, min = 0.001, max = 0.5)
                   )
            )
          ),

          fluidRow(
            column(12,
                   div(style = "text-align: center; margin: 30px 0;",
                       actionButton("run_fisher_bound", "Calculate Fisher Bound",
                                    class = "btn-primary btn-lg",
                                    icon = icon("border-all"),
                                    style = "padding: 15px 50px; font-size: 1.3em;")
                   ),
                   div(id = "fisher_bound_status_message", style = "text-align: center; margin-top: 20px;")
            )
          ),

          fluidRow(
            column(12,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08); margin-bottom: 20px;",
                       h4(icon("table"), " Fisher Bound Results",
                          style = "color: #667eea; font-weight: 600; margin-bottom: 20px;"),
                       fluidRow(
                         column(6,
                                div(style = "background: #f8f9fc; padding: 15px; border-radius: 6px; border-left: 4px solid #667eea; margin-bottom: 15px;",
                                    h5("Contingency Table (M)", style = "color: #667eea; margin-bottom: 10px;"),
                                    verbatimTextOutput("fisherBoundM")
                                )
                         ),
                         column(6,
                                div(style = "background: #f8f9fc; padding: 15px; border-radius: 6px; border-left: 4px solid #764ba2; margin-bottom: 15px;",
                                    h5("P-value at Boundary (p)", style = "color: #764ba2; margin-bottom: 10px;"),
                                    verbatimTextOutput("fisherBoundP")
                                )
                         )
                       ),
                       fluidRow(
                         column(6,
                                div(style = "background: #f8f9fc; padding: 15px; border-radius: 6px; border-left: 4px solid #11998e; margin-bottom: 15px;",
                                    h5("Responders for Control Arm (rc)", style = "color: #11998e; margin-bottom: 10px;"),
                                    verbatimTextOutput("fisherBoundRc")
                                )
                         ),
                         column(6,
                                div(style = "background: #f8f9fc; padding: 15px; border-radius: 6px; border-left: 4px solid #f093fb; margin-bottom: 15px;",
                                    h5("Sample Size in Control Arm (nc)", style = "color: #f093fb; margin-bottom: 10px;"),
                                    verbatimTextOutput("fisherBoundNc")
                                )
                         )
                       ),
                       fluidRow(
                         column(6,
                                div(style = "background: #f8f9fc; padding: 15px; border-radius: 6px; border-left: 4px solid #667eea; margin-bottom: 15px;",
                                    h5("Response Rate for Experimental Arm (rt)", style = "color: #667eea; margin-bottom: 10px;"),
                                    verbatimTextOutput("fisherBoundRt")
                                )
                         ),
                         column(6,
                                div(style = "background: #f8f9fc; padding: 15px; border-radius: 6px; border-left: 4px solid #764ba2; margin-bottom: 15px;",
                                    h5("Sample Size in Experimental Arm (nt)", style = "color: #764ba2; margin-bottom: 10px;"),
                                    verbatimTextOutput("fisherBoundNt")
                                )
                         )
                       ),
                       fluidRow(
                         column(12,
                                div(style = "background: #f8f9fc; padding: 15px; border-radius: 6px; border-left: 4px solid #11998e;",
                                    h5("Minimum Detectable Difference (delta)", style = "color: #11998e; margin-bottom: 10px;"),
                                    verbatimTextOutput("fisherBoundDelta")
                                )
                         )
                       )
                   )
            )
          ),

          fluidRow(
            column(12,
                   div(style = "background: white; padding: 25px; border-radius: 8px; box-shadow: 0 4px 15px rgba(0,0,0,0.08);",
                       h4(icon("lightbulb"), " Interpretation",
                          style = "color: #667eea; font-weight: 600; margin-bottom: 15px;"),
                       htmlOutput("fisherBoundConclusion"),
                       htmlOutput("fisherBoundStatisticalConclusion")
                   )
            )
          )
        )
      )
    )
  )
)

