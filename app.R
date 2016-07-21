library(shiny)
library(ggplot2)

ui <- fluidPage(sidebarLayout(
   sidebarPanel(
      
      # Variables
      
      tags$h3("Input Variables"),
      sliderInput(
         "n",
         "Manning's n",
         min = 0.000,
         max = 0.300,
         value = c(0.03, 0.04)
      ),
      sliderInput(
         "WidthT",
         "Top Width (m)",
         min = 1,
         max = 1000,
         value = c(10, 20)
      ),
      sliderInput(
         "WidthB",
         "Bottom Width (m)",
         min = 1,
         max = 1000,
         value = c(5, 10)
      ),
      sliderInput(
         "Depth",
         "Depth (m)",
         min = 0.2,
         max = 100,
         value = c(1, 2)
      ),
      sliderInput(
         "SlopeB",
         "Bed Slope (m/m)",
         min = 0.0001,
         max = 0.25,
         value = c(0.002, 0.005)
      ),
      
      downloadButton("mycsv", "Download csv file")
   ),
   mainPanel(
      tags$h2("Monte Carlo Analysis of Manning's Equation"),
      
      #tableOutput("table"),
      "An app demonstrating impacts of input value uncertainty on results of Manning's equation for open channel flow using Monte Carlo analysis.  See notes at bottom of screen for details.",
      
      plotOutput("box"),
      verbatimTextOutput("stats"),

      "Programmed by John Yagecic, P.E. in US customary units (JYagecic@gmail.com)",
      tags$br(),
      "Adapted to metric units by Mike Spencer (http://mikerspencer.com)",
      tags$br(),
      tags$br(),
      tags$a(href = "https://en.wikipedia.org/wiki/Manning_formula", "More about Manning's Equation."),
      tags$br(),
      tags$a(href = "https://en.wikipedia.org/wiki/Monte_Carlo_method", "More about Monte Carlo method"),
      plotOutput("fourpanel"),
      #verbatimTextOutput("ManningDF"),
      "All distributions are uniform distributions with minimum and maximum values set by the slider bars.",
      "Input variable vectors have a length of 10,000 corresponding to 10,000 Monte Carlo iterations.",
      "Computation assumes a symmetrical trapezoidal channel.",
      "Bottom width should be less than top width for reasonable results.",
      tags$br(),
      "If you use this product or the underlying code in any professional or academic product, please consider ",
      "using a citation such as:",
      tags$br(),
      tags$br(),
      "Yagecic, John, July 2016.  Monte Carlo Analysis of Manning's Equation: a web app demonstrating impacts of input value uncertainty on results of Manning's equation for open channel flow using Monte Carlo analysis.",
      tags$br(),
      tags$br(),
      tags$a(href = "https://github.com/mikerspencer", "Get the script in metric units"),
      tags$br(),
      tags$a(href = "https://github.com/JohnYagecic", "Get the script in US customary units")
   )
))

server <- function(input, output) {
   n <- reactive({
      runif(10000, input$n[1], input$n[2])
   })
   WidthT <-
      reactive({
         runif(10000, input$WidthT[1], input$WidthT[2])
      })
   WidthB <-
      reactive({
         runif(10000, input$WidthB[1], input$WidthB[2])
      })
   Depth <-
      reactive({
         runif(10000, input$Depth[1], input$Depth[2])
      })
   SlopeB <-
      reactive({
         runif(10000, input$SlopeB[1], input$SlopeB[2])
      })
   Adj <- reactive({
      (WidthT() - WidthB()) / 2
   })
   Opp <- reactive({
      Depth()
   })
   Hyp <- reactive({
      sqrt(Adj() ^ 2 + Opp() ^ 2)
   })
   AreaSide <- reactive({
      0.5 * Adj() * Opp()
   })
   AreaTotal <- reactive({
      (Depth() * WidthB()) + 2 * AreaSide()
   })
   WettedPerim <- reactive({
      WidthB() + 2 * Hyp()
   })
   RadiusHyd <- reactive({
      AreaTotal() / WettedPerim()
   })
   Vel <-
      reactive({
         (1 / n()) * ((RadiusHyd()) ^ (2 / 3)) * (SlopeB() ^ 0.5)
      })
   Q <- reactive({
      Vel() * AreaTotal()
   })
   
   ManningDF <- reactive({
      data.frame(
         n = n(),
         TopWidth = WidthT(),
         BottomWidth = WidthB(),
         Depth = Depth(),
         BedSlope = SlopeB(),
         Area = AreaTotal(),
         WettedPerimeter = WettedPerim(),
         HydraulicRadius = RadiusHyd(),
         Velocity = Vel(),
         Discharge = Q()
      )
   })
   
   
   
   output$box <- renderPlot({
      par(mfrow = c(1, 2))
      hist(Q(),
           breaks = 20,
           main = "Histogram of Q (CUMECS)",
           xlab = "Q (CUMECS)")
      boxplot(Q(), ylab = "Q (CUMECS)", main = "Boxplot of Q (CUMECS)")
   })
   output$stats <- renderPrint({
      summary(Q())
   })
   
   
   output$fourpanel <- renderPlot({
      par(mfrow = c(2, 2))
      hist(
         AreaTotal(),
         main = "Histogram of Area",
         xlab = "Area (sq. m)",
         breaks = 20,
         col = "blue"
      )
      hist(
         Vel(),
         main = "Histogram of Velocity",
         xlab = "Velocity (m/s)",
         breaks = 20,
         col = "red"
      )
      hist(
         WettedPerim(),
         main = "Histogram of Wetted Perminter",
         xlab = "Wetted Perimeter (m)",
         breaks = 20,
         col = "green"
      )
      hist(
         RadiusHyd(),
         main = "Histogram of Hydraulic Radius",
         xlab = "Hydraulic Radius",
         breaks = 20,
         col = "purple"
      )
   })
   
   #output$ManningDF <-renderPrint({head(ManningDF())})
   
   output$mycsv <- downloadHandler(
      filename = c('ManningMCdata.csv'),
      content = function(file) {
         setwd(tempdir())
         write.csv(ManningDF(), file)
      }
   )
   
   
}

shinyApp(ui = ui, server = server)