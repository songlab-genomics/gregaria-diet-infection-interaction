#!/usr/bin/env node
// Build the manuscript supplementary workbook from validated JSON tables.

import fs from "node:fs/promises";
import path from "node:path";
import {
  FileBlob,
  SpreadsheetFile,
  Workbook,
} from "@oai/artifact-tool";

function argument(name) {
  const index = process.argv.indexOf(name);
  if (index < 0 || index + 1 >= process.argv.length) {
    throw new Error(`Missing required argument: ${name}`);
  }
  return process.argv[index + 1];
}

const projectDir = path.resolve(argument("--project"));
const jsonPath = path.resolve(argument("--json"));
const outputPath = path.resolve(argument("--output"));
const previewDir = path.resolve(argument("--preview-dir"));

await fs.mkdir(path.dirname(outputPath), { recursive: true });
await fs.mkdir(previewDir, { recursive: true });

const payload = JSON.parse(await fs.readFile(jsonPath, "utf8"));

// Render the manual-curation workbook first so its structure is checked before use.
const candidateWorkbookPath = path.join(
  projectDir,
  "data/external/curated_target_genes/2026-08-09_Immune_Protein_Anabolism_DEGs_v24.xlsx",
);
const sourceWorkbook = await SpreadsheetFile.importXlsx(
  await FileBlob.load(candidateWorkbookPath),
);
for (const [sheetName, range] of [
  ["Immune", "A1:V25"],
  ["Protein-Anabolism", "A1:V25"],
]) {
  const preview = await sourceWorkbook.render({
    sheetName,
    range,
    scale: 1,
    format: "png",
  });
  await fs.writeFile(
    path.join(previewDir, `source_${sheetName.replaceAll("-", "_")}.png`),
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const workbook = Workbook.create();
const headerFill = "#263238";
const headerFont = "#FFFFFF";
const borderColor = "#D8DEE3";
const titleFill = "#E8F1F2";
const noteFill = "#F5F7F8";

function excelColumnName(index) {
  let value = index + 1;
  let name = "";
  while (value > 0) {
    const remainder = (value - 1) % 26;
    name = String.fromCharCode(65 + remainder) + name;
    value = Math.floor((value - 1) / 26);
  }
  return name;
}

function safeTableName(sheetName) {
  return `${sheetName.replaceAll(/[^A-Za-z0-9]/g, "")}Table`;
}

function columnWidth(header) {
  const text = header.toLowerCase();
  if (text.includes("description") || text.includes("annotation") || text.includes("reason")) return 330;
  if (text.includes("geneid") || text.includes("significant_contrasts")) return 300;
  if (text.includes("pathway") || text.includes("pfam") || text === "gos") return 280;
  if (text.includes("source") || text.includes("notes")) return 300;
  if (text.includes("comparison") || text.includes("condition")) return 220;
  if (text.includes("theme") || text.includes("category") || text.includes("evidence")) return 220;
  if (text.includes("status")) return 250;
  if (text.includes("gene_id") || text === "id" || text.includes("term_id")) return 135;
  if (text.includes("label") || text.includes("sample")) return 120;
  return 115;
}

function applyColumnFormatting(sheet, headers, rowCount) {
  headers.forEach((header, columnIndex) => {
    const column = sheet.getRangeByIndexes(0, columnIndex, rowCount + 1, 1);
    column.format.columnWidthPx = columnWidth(header);
    const text = header.toLowerCase();

    if (
      text.includes("description") || text.includes("annotation") ||
      text.includes("reason") || text.includes("geneid") ||
      text.includes("pathway") || text.includes("pfam") || text === "gos" ||
      text.includes("source") || text.includes("notes") ||
      text.includes("significant_contrasts")
    ) {
      column.format.wrapText = true;
      column.format.verticalAlignment = "top";
    }
    if (text.includes("padj") || text.includes("p.adjust") || text.includes("pvalue") || text === "qvalue") {
      column.format.numberFormat = "0.00E+00";
    } else if (text.includes("log2foldchange") || text.includes("log2fc") || text.includes("foldenrichment") || text.includes("zscore")) {
      column.format.numberFormat = "0.00";
    } else if (text.endsWith("_pct") || text.includes("mapping_pct")) {
      column.format.numberFormat = "0.0000";
    } else if (text.includes("reads") || text === "count" || text.startsWith("n_")) {
      column.format.numberFormat = "#,##0";
    }
  });
}

function addConditionalFormatting(sheet, headers, rowCount) {
  const headerIndex = Object.fromEntries(headers.map((header, index) => [header, index]));
  for (const logColumn of ["log2FoldChange", "maximum_significant_abs_log2fc"]) {
    if (headerIndex[logColumn] === undefined) continue;
    const letter = excelColumnName(headerIndex[logColumn]);
    const range = sheet.getRange(`${letter}2:${letter}${rowCount + 1}`);
    range.conditionalFormats.add("cellIs", {
      operator: "greaterThan",
      formula: 0,
      format: { font: { color: "#B91C1C", bold: true } },
    });
    if (logColumn === "log2FoldChange") {
      range.conditionalFormats.add("cellIs", {
        operator: "lessThan",
        formula: 0,
        format: { font: { color: "#1D4ED8", bold: true } },
      });
    }
  }

  if (headerIndex.evidence_tier !== undefined) {
    const letter = excelColumnName(headerIndex.evidence_tier);
    const range = sheet.getRange(`${letter}2:${letter}${rowCount + 1}`);
    range.conditionalFormats.add("containsText", {
      text: "Manual curation candidate",
      format: { fill: "#C9D8E8", font: { bold: true, color: "#17324D" } },
    });
    range.conditionalFormats.add("containsText", {
      text: "eggNOG-supported candidate",
      format: { fill: "#EFE8F5", font: { color: "#6B4C7A" } },
    });
  }

  if (headerIndex.analysis_included !== undefined) {
    const letter = excelColumnName(headerIndex.analysis_included);
    const range = sheet.getRange(`${letter}2:${letter}${rowCount + 1}`);
    range.conditionalFormats.add("cellIs", {
      operator: "equal",
      formula: "FALSE",
      format: { fill: "#FCE8E6", font: { bold: true, color: "#A61B1B" } },
    });
  }
}

function addDataSheet(sheetName, tablePayload, style = "TableStyleMedium2") {
  const sheet = workbook.worksheets.add(sheetName);
  sheet.showGridLines = false;
  const headers = tablePayload.columns;
  const rows = tablePayload.rows;
  const matrix = [headers, ...rows];
  const range = sheet.getRangeByIndexes(0, 0, matrix.length, headers.length);
  range.values = matrix;
  range.format.font = { size: 10, color: "#263238" };
  range.format.rowHeight = 20;
  range.format.borders = {
    insideHorizontal: { style: "thin", color: borderColor },
  };
  const header = sheet.getRangeByIndexes(0, 0, 1, headers.length);
  header.format = {
    fill: headerFill,
    font: { bold: true, color: headerFont, size: 10 },
    wrapText: true,
    verticalAlignment: "center",
    rowHeight: 36,
  };
  applyColumnFormatting(sheet, headers, rows.length);
  addConditionalFormatting(sheet, headers, rows.length);
  const finalCell = `${excelColumnName(headers.length - 1)}${rows.length + 1}`;
  const table = sheet.tables.add(`A1:${finalCell}`, true, safeTableName(sheetName));
  table.style = style;
  table.showFilterButton = true;
  sheet.freezePanes.freezeRows(1);
  sheet.freezePanes.freezeColumns(Math.min(2, headers.length));
  return sheet;
}

const readme = workbook.worksheets.add("README");
readme.showGridLines = false;
readme.getRange("A1:H2").merge();
readme.getRange("A1").values = [[payload.metadata.title]];
readme.getRange("A1:H2").format = {
  fill: titleFill,
  font: { bold: true, color: "#17324D", size: 20 },
  verticalAlignment: "center",
};
readme.getRange("A4:B9").values = [
  ["Analysis", payload.metadata.analysis],
  ["Cohort", payload.metadata.cohort],
  ["Gene scope", payload.metadata.gene_scope],
  ["Candidate evidence rule", payload.metadata.candidate_rule],
  ["rRNA IDs in DEG catalogue", payload.metadata.rrna_overlap_in_deg_catalogue],
  ["Workbook use", "Filter any data sheet by gene, contrast, category, treatment, or adjusted P value."],
];
readme.getRange("A4:A9").format = {
  fill: headerFill,
  font: { bold: true, color: headerFont },
};
readme.getRange("B4:B9").format = {
  fill: noteFill,
  wrapText: true,
  verticalAlignment: "top",
};
readme.getRange("A4:B9").format.borders = {
  preset: "all",
  style: "thin",
  color: borderColor,
};
readme.getRange("A4:A9").format.columnWidthPx = 190;
readme.getRange("B4:B9").format.columnWidthPx = 680;
readme.getRange("A11:H11").merge();
readme.getRange("A11").values = [[
  "Dark evidence cells identify manual curation candidates. Lighter evidence cells identify eggNOG-supported candidates; the evidence and matched rule remain explicit columns.",
]];
readme.getRange("A11:H11").format = {
  fill: "#FFF7E6",
  font: { italic: true, color: "#5F4B1C" },
  wrapText: true,
  rowHeight: 42,
};
readme.freezePanes.freezeRows(2);

addDataSheet("Sources", payload.tables.Sources, "TableStyleLight9");
addDataSheet("Samples_mapping", payload.tables.Samples_mapping, "TableStyleMedium4");
addDataSheet("DEG_transcript_exon", payload.tables.DEG_transcript_exon, "TableStyleMedium2");
addDataSheet("GO_enrichment", payload.tables.GO_enrichment, "TableStyleMedium4");
addDataSheet("KEGG_enrichment", payload.tables.KEGG_enrichment, "TableStyleMedium9");
addDataSheet("Physiology_candidates", payload.tables.Physiology_candidates, "TableStyleMedium2");
addDataSheet("Candidate_DEG_results", payload.tables.Candidate_DEG_results, "TableStyleMedium2");

const summaryInspect = await workbook.inspect({
  kind: "workbook,sheet,table",
  maxChars: 9000,
  tableMaxRows: 4,
  tableMaxCols: 8,
});
console.log(summaryInspect.ndjson);

const errorInspect = await workbook.inspect({
  kind: "match",
  searchTerm: "#REF!|#DIV/0!|#VALUE!|#NAME\\?|#N/A",
  options: { useRegex: true, maxResults: 200 },
  summary: "final formula error scan",
});
console.log(errorInspect.ndjson);

const previewRanges = {
  README: "A1:H11",
  Sources: "A1:D8",
  Samples_mapping: "A1:R16",
  DEG_transcript_exon: "A1:W16",
  GO_enrichment: "A1:R16",
  KEGG_enrichment: "A1:S16",
  Physiology_candidates: "A1:T18",
  Candidate_DEG_results: "A1:T18",
};
for (const [sheetName, range] of Object.entries(previewRanges)) {
  const preview = await workbook.render({
    sheetName,
    range,
    scale: 1,
    format: "png",
  });
  await fs.writeFile(
    path.join(previewDir, `${sheetName}.png`),
    new Uint8Array(await preview.arrayBuffer()),
  );
}

const output = await SpreadsheetFile.exportXlsx(workbook);
await output.save(outputPath);
console.log(`Saved ${outputPath}`);
