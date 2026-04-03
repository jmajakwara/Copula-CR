library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Create the publication year trend data
year_trend <- data.frame(
  Year = c(1994, 1995, 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 
           2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011, 2012, 2013,
           2014, 2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022, 2023,
           2024, 2025, 2026),
  Articles = c(2, 2, 0, 0, 0, 0, 0, 1, 1, 4,
               2, 3, 4, 4, 3, 7, 4, 8, 9, 5,
               11, 6, 13, 11, 17, 12, 5, 18, 18, 10,
               11, 5, 2)
)
# Calculate cumulative
year_trend$Cumulative <- cumsum(year_trend$Articles)

####  CUMULATIVE PUBLICATIONS LINE CHART
# ============================================================

p <- ggplot(year_trend, aes(x = Year, y = Cumulative)) +
  # Line for cumulative publications
  geom_line(color = "darkgreen", linewidth = 1.5) +
  geom_point(color = "darkgreen", size = 3, alpha = 0.7) +
  
  # Add area fill
  geom_area(fill = "lightgreen", alpha = 0.4) +
  
  # Labels and titles
  labs(
  #  title = "Cumulative Growth of Research in Copula-Based CR Methods",
  #  subtitle = "Total Publications from 1994 to 27 February 2026",
    x = "Publication Year",
    y = "Cumulative Number of Articles"
 #   caption = "Total of 186 articles as of 27 February 2026"
  ) +
  
  # Scale x-axis
  scale_x_continuous(breaks = seq(1994, 2026, by = 4), 
                     minor_breaks = seq(1994, 2026, by = 2)) +
  
  # Scale y-axis
  scale_y_continuous(limits = c(0, 220), breaks = seq(0, 220, by = 25)) +
  
  # Theme customization
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    plot.caption = element_text(size = 9, color = "gray50", hjust = 1),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_line(color = "gray90"),
    panel.grid.major = element_line(color = "gray80")
  )

print(p)
