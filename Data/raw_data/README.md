# Raw data

The field data are not distributed with this repository. They are available from
the corresponding author on reasonable request.

To run the analysis, place the spreadsheet in this folder as

    Data/raw_data/Dados_Percevejos_Antiteuchus_Monalonion.xlsx

The script `R_Scripts/00_Data_Preparation.R` reads the sheet `Planilha1`, which
has one row per sampled tree and the columns

| Column      | Content                                            |
|-------------|----------------------------------------------------|
| ID          | tree identifier (each tree was sampled only once)  |
| Date        | sampling date                                      |
| Month       | sampling month                                     |
| Time        | days since the first sampling date                 |
| Altitude    | elevation class: High, Medium or Low               |
| Farm        | farm identifier                                    |
| Condition   | shading condition: Shaded or Unshaded              |
| Monalonion  | number of *Monalonion* individuals on the tree     |
| Antiteuchus | number of *Antiteuchus* individuals on the tree    |
| Temperature | air temperature                                    |
| Humidity    | relative humidity                                  |
