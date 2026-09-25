import { computeSheet, readNumber } from "./calc.mjs";
import { ms2WorkedExample, t4ReportedWater } from "./examples.mjs";

const STORAGE_KEY = "amb-field-sheet-v1";

const SCALAR_FIELDS = [
  ["herb", "herb"],
  ["virus", "virus"],
  ["unitSystem", "unitSystem"],
  ["operator", "operator"],
  ["date", "date"],
  ["hours", "hours"],
  ["notes", "notes"],
  ["stockConcentration", "stockConcentration"],
  ["stockVolumeMl", "stockVolumeMl"],
  ["diluentVolumeMl", "diluentVolumeMl"],
  ["cultureVolumeMl", "cultureVolumeMl"],
  ["blankDecayPercent", "blankDecayPercent"],
  ["blankInitial", "blankInitial"],
  ["blankFinal", "blankFinal"],
  ["waterConcentration", "waterConcentrationPerMl"],
  ["genomeLength", "genomeLengthBp"],
  ["molecularWeight", "molecularWeight"],
  ["targetParticles", "targetParticles"],
  ["rootParticles", "rootParticles"],
  ["waterParticles", "waterParticles"],
  ["mixtureA", "mixtureA"],
  ["mixtureExtract", "mixtureExtractMl"],
  ["mixtureSum", "mixtureSumFwtG"],
  ["mixtureHerb", "mixtureWholeHerbG"],
];

const PART_FIELDS = [
  ["name", "name"],
  ["mass", "massG"],
  ["plaque", "plaqueCount"],
  ["dilution", "dilutionFactor"],
  ["inoculum", "inoculumMl"],
  ["extract", "extractMl"],
  ["infective", "infectiveTotal"],
  ["inactivated", "inactivatedTotal"],
  ["fraction", "inactivatedFraction"],
  ["dna", "dnaNgPerUl"],
];

let currentId = null;
let latestInput = null;
let latestResult = null;
let lastErrorSignature = "";
let clearArmed = false;

function $(id) {
  return document.getElementById(id);
}

function today() {
  return new Date().toISOString().slice(0, 10);
}

function emptyInput() {
  return {
    herb: "",
    virus: "ms2",
    unitSystem: "PFU",
    operator: "",
    date: today(),
    hours: "",
    notes: "",
    stockConcentration: "",
    stockVolumeMl: "",
    diluentVolumeMl: "",
    cultureVolumeMl: "",
    blankMode: "percent",
    blankDecayPercent: "",
    blankInitial: "",
    blankFinal: "",
    waterConcentrationPerMl: "",
    genomeLengthBp: 3569,
    molecularWeight: 330,
    targetParticles: "",
    rootParticles: "",
    waterParticles: "",
    mixtureA: "",
    mixtureExtractMl: "",
    mixtureSumFwtG: "",
    mixtureWholeHerbG: "",
    parts: [{ name: "root" }, { name: "leaf" }, { name: "shoot" }],
  };
}

function inputNumber(value) {
  if (value === null || value === undefined || value === "") return "";
  const parsed = typeof value === "number" ? value : Number(value);
  if (!Number.isFinite(parsed)) return String(value);
  const abs = Math.abs(parsed);
  if (abs !== 0 && (abs >= 1e6 || abs < 1e-4)) return parsed.toExponential(8);
  return String(parsed);
}

function formatCount(value) {
  if (!Number.isFinite(value)) return "n/a";
  const abs = Math.abs(value);
  if (abs !== 0 && (abs < 0.01 || abs >= 10000)) return value.toExponential(3);
  return new Intl.NumberFormat("en-US", { maximumFractionDigits: 3 }).format(value);
}

function formatPercent(value) {
  if (!Number.isFinite(value)) return "n/a";
  return `${value.toFixed(1)}%`;
}

function withUnit(value, unit) {
  const text = formatCount(value);
  return text === "n/a" ? text : `${text} ${unit}`;
}

function virusLabel(virus) {
  if (virus === "t4") return "T4";
  if (virus === "other") return "Other virus";
  return "MS2";
}

function sheetTitle(input) {
  const bits = [String(input.herb || "").trim(), virusLabel(input.virus)];
  if (readNumber(input.hours) !== null) bits.push(`${input.hours} h`);
  const title = bits.filter((bit) => bit).join(" / ");
  return title || "Untitled record";
}

function setStatus(message) {
  $("status").textContent = message;
}

function applyVirusPreset() {
  const virus = $("virus").value;
  const length = $("genomeLength");
  const mass = $("molecularWeight");
  const presets = {
    ms2: ["3569", "330"],
    t4: ["169000", "660"],
  };
  if (!presets[virus]) {
    length.readOnly = false;
    mass.readOnly = false;
    return;
  }
  [length.value, mass.value] = presets[virus];
  length.readOnly = true;
  mass.readOnly = true;
}

function applyBlankMode() {
  const pair = document.querySelector('input[name="blankMode"]:checked')?.value === "pair";
  $("blank-percent-fields").hidden = pair;
  $("blank-pair-fields").hidden = !pair;
}

function syncLegend(part) {
  const name = part.querySelector('[data-field="name"]').value.trim();
  part.querySelector(".legend-name").textContent = name || "Plant part";
}

function renumberParts() {
  const parts = document.querySelectorAll("[data-part]");
  parts.forEach((part, index) => {
    part.dataset.index = String(index);
    part.querySelectorAll("[data-error-key]").forEach((field) => {
      field.dataset.errorId = `part-${index}-${field.dataset.errorKey}`;
    });
    const remove = part.querySelector("[data-remove]");
    remove.disabled = parts.length === 1;
    syncLegend(part);
  });
}

function createPart(values = {}) {
  const node = $("part-template").content.firstElementChild.cloneNode(true);
  for (const [field, key] of PART_FIELDS) {
    const input = node.querySelector(`[data-field="${field}"]`);
    const value = values[key];
    input.value = typeof value === "number" ? inputNumber(value) : value == null ? "" : String(value);
  }
  node.querySelector('[data-field="include"]').checked = Boolean(values.includeHdp);
  syncLegend(node);
  return node;
}

function setParts(parts) {
  const host = $("parts");
  host.replaceChildren();
  const rows = parts && parts.length ? parts : [{ name: "root" }];
  for (const part of rows) host.append(createPart(part));
  renumberParts();
}

function readPart(node) {
  const value = (field) => node.querySelector(`[data-field="${field}"]`).value;
  return {
    name: value("name"),
    massG: value("mass"),
    plaqueCount: value("plaque"),
    dilutionFactor: value("dilution"),
    inoculumMl: value("inoculum"),
    extractMl: value("extract"),
    infectiveTotal: value("infective"),
    inactivatedTotal: value("inactivated"),
    inactivatedFraction: value("fraction"),
    dnaNgPerUl: value("dna"),
    includeHdp: node.querySelector('[data-field="include"]').checked,
  };
}

function readForm() {
  const input = { parts: [] };
  for (const [id, key] of SCALAR_FIELDS) input[key] = $(id).value;
  input.blankMode = document.querySelector('input[name="blankMode"]:checked')?.value || "percent";
  input.parts = [...document.querySelectorAll("[data-part]")].map(readPart);
  return input;
}

function fillForm(input) {
  for (const [id, key] of SCALAR_FIELDS) {
    if (!(key in input)) continue;
    const value = input[key];
    $(id).value = typeof value === "number" ? inputNumber(value) : value == null ? "" : String(value);
  }
  const mode = input.blankMode === "pair" ? "pair" : "percent";
  const radio = document.querySelector(`input[name="blankMode"][value="${mode}"]`);
  if (radio) radio.checked = true;
  applyVirusPreset();
  if ($("virus").value === "other") {
    if ("genomeLengthBp" in input) $("genomeLength").value = inputNumber(input.genomeLengthBp);
    if ("molecularWeight" in input) $("molecularWeight").value = inputNumber(input.molecularWeight);
  }
  applyBlankMode();
  setParts(input.parts || [{ name: "root" }]);
  recompute();
}

function markErrors(result) {
  document.querySelectorAll("[data-error-id]").forEach((field) => field.removeAttribute("aria-invalid"));
  const ids = new Set(result.errors.map((error) => error.id));
  document.querySelectorAll("[data-error-id]").forEach((field) => {
    if (ids.has(field.dataset.errorId)) field.setAttribute("aria-invalid", "true");
  });
  const signature = result.errors.map((error) => `${error.id}:${error.message}`).join("|");
  const list = $("result-errors");
  if (signature === lastErrorSignature) return;
  lastErrorSignature = signature;
  list.replaceChildren();
  for (const error of result.errors) {
    const item = document.createElement("li");
    item.textContent = displayError(error.message);
    list.append(item);
  }
}

function addLine(list, term, value, id) {
  const row = document.createElement("div");
  const dt = document.createElement("dt");
  dt.textContent = term;
  const dd = document.createElement("dd");
  dd.textContent = value;
  if (id) dd.id = id;
  row.append(dt, dd);
  list.append(row);
}

function joinWords(items) {
  if (items.length <= 1) return items[0] || "";
  if (items.length === 2) return `${items[0]} and ${items[1]}`;
  return `${items.slice(0, -1).join(", ")} and ${items.at(-1)}`;
}

const PYTHON_MESSAGES = {
  "initial_infectivity must be positive": "The initial value must be greater than zero.",
  "dilution_factor and virus_diluted_volume_ml must be positive":
    "The dilution and the inoculum volume must be greater than zero.",
  "fresh_weight_g must be positive": "Fresh weight must be greater than zero.",
  "pi_by_part and mass_by_part must share at least one key":
    "Combined infectivity requires a matching infectivity and fresh weight.",
  "sum of masses must be positive": "The sum of fresh weights must be greater than zero.",
  "sum_fwt_parts_g must be positive": "The sum of fresh weights must be greater than zero.",
  "initial_particles must be positive": "The initial concentration must be greater than zero.",
  "viral_particles_root must be positive": "The root particle count must be greater than zero.",
  "viral_particles_water must be positive": "The water particle count must be greater than zero.",
  "target_length_bp and molecular_weight_eq9 must be positive":
    "Genome length and the mass term must be greater than zero.",
  "n0_eff_adj must be positive": "The adjusted initial amount must be greater than zero.",
  "total volume must be positive": "Total volume must be greater than zero.",
};

function displayError(message) {
  return PYTHON_MESSAGES[message] || message;
}

function resultLead(result) {
  if (result.errors.length) {
    return "Correct the marked fields. Results that can be calculated are still shown.";
  }
  if (!result.balance) {
    const names = {
      spike: "the initial spike",
      blank: "the blank decay",
      water: "the water concentration",
    };
    const missing = result.missing.map((key) => names[key]).filter(Boolean);
    if (missing.length) return `Enter ${joinWords(missing)} to close the mass balance.`;
    return "A harvested part does not have an infective total.";
  }
  const residual = result.balance.residualFraction;
  let lead = "The water and the plant parts account for the adjusted initial amount.";
  if (residual > 0) {
    lead = "Part of the adjusted initial amount is not recovered in the water or the plant parts.";
  }
  if (residual < 0) lead = "The water and the plant parts exceed the adjusted initial amount.";
  if (result.hdpIncluded) {
    return `${lead} Highly degraded particle totals are included in this closure. Include them only when they use the same unit as the initial spike.`;
  }
  const showedHdp = result.parts.some((part) => part.hdpPerMl !== null);
  if (showedHdp) {
    return `${lead} Highly degraded particles are calculated and excluded from the closure.`;
  }
  return lead;
}

function renderResults(input, result) {
  const unit = input.unitSystem || "PFU";
  markErrors(result);
  $("result-lead").textContent = resultLead(result);
  $("spike-hint").textContent =
    result.c0 === null
      ? "The initial concentration is shown after the stock concentration and both volumes are entered."
      : `Initial concentration (C0) is ${withUnit(result.c0, unit)}. Culture volume is ${formatCount(result.cultureVolumeMl)} mL.`;
  $("blank-hint").textContent =
    result.blankDecayPercent === null
      ? "Natural decay from the herb-free blank is subtracted before closure."
      : `Blank decay is ${formatPercent(result.blankDecayPercent)}. This value is subtracted before closure.`;
  $("water-hint").textContent =
    result.pePercent === null
      ? "Phytoremediation efficiency is shown after the initial and final water concentrations are entered."
      : `Phytoremediation efficiency is ${formatPercent(result.pePercent)}.`;

  const numbers = $("result-numbers");
  numbers.replaceChildren();
  addLine(numbers, "Initial concentration (C0)", withUnit(result.c0, `${unit}/mL`), "panel-c0");
  addLine(numbers, "Blank decay", formatPercent(result.blankDecayPercent), "panel-blank");
  addLine(numbers, "Phytoremediation efficiency", formatPercent(result.pePercent), "panel-pe");
  if (result.balance) {
    addLine(numbers, "Adjusted initial amount", withUnit(result.balance.n0EffAdj, unit), "panel-n0");
    addLine(numbers, "Water total", withUnit(result.balance.nWater, unit), "panel-water");
    addLine(numbers, "Plant total", withUnit(result.balance.nPlantTotal, unit), "panel-plant");
    addLine(numbers, "Residual", withUnit(result.balance.nResidual, unit), "panel-residual-count");
    addLine(numbers, "Closure", formatPercent(result.balance.closureFraction * 100), "panel-closure");
    addLine(numbers, "Residual fraction", formatPercent(result.balance.residualFraction * 100), "panel-residual");
  } else {
    addLine(numbers, "Closure", "n/a", "panel-closure");
    addLine(numbers, "Residual fraction", "n/a", "panel-residual");
  }
  if (result.tf !== null) addLine(numbers, "Translocation factor", formatCount(result.tf), "panel-tf");
  if (result.bcf !== null) addLine(numbers, "Bioconcentration factor", formatCount(result.bcf), "panel-bcf");
  const weighted = result.parts.filter((part) => part.used && part.pi !== null && part.massG > 0);
  if (weighted.length >= 2 && result.combinedPi !== null) {
    addLine(numbers, "Combined infectivity (Eq. 4)", withUnit(result.combinedPi, `${unit}/g`), "panel-eq4");
  }
  if (result.mixturePi !== null) {
    addLine(numbers, "Mixture infectivity (Eq. 5)", withUnit(result.mixturePi, `${unit}/g`), "panel-eq5");
  }

  const parts = $("result-parts");
  parts.replaceChildren();
  for (const part of result.parts.filter((row) => row.used)) {
    const block = document.createElement("div");
    block.className = "part-result";
    const title = document.createElement("p");
    title.className = "sheet-kicker";
    title.textContent = part.name;
    block.append(title);
    const lines = [];
    if (part.a !== null) lines.push(`Infectivity (A) ${withUnit(part.a, `${unit}/mL`)}`);
    if (part.pi !== null) lines.push(`Infectivity per gram (PI) ${withUnit(part.pi, `${unit}/g`)}`);
    if (part.infective !== null) lines.push(`Infective total ${withUnit(part.infective, unit)}`);
    if (part.inactivated !== null) lines.push(`Inactivated total ${withUnit(part.inactivated, unit)}`);
    if (part.hdpPerMl !== null) {
      lines.push(`Highly degraded particles ${withUnit(part.hdpPerMl, "Vp/mL")}`);
    }
    if (part.hdpTotal !== null) {
      lines.push(
        part.includedHdp
          ? `Highly degraded particle total ${withUnit(part.hdpTotal, "Vp")} included in the closure`
          : `Highly degraded particle total ${withUnit(part.hdpTotal, "Vp")} excluded from the closure`,
      );
    }
    for (const line of lines) {
      const row = document.createElement("p");
      row.textContent = line;
      block.append(row);
    }
    parts.append(block);
  }

  $("notebook").textContent = notebookText(input, result, unit);
  $("bar-closure").textContent = result.balance ? formatPercent(result.balance.closureFraction * 100) : "n/a";
  $("bar-residual").textContent = result.balance ? formatPercent(result.balance.residualFraction * 100) : "n/a";
  $("bar-pe").textContent = formatPercent(result.pePercent);
}

function notebookText(input, result, unit) {
  const lines = [sheetTitle(input)];
  if (input.date) lines.push(input.date);
  if (input.operator) lines.push(input.operator);
  lines.push(`Measurement unit ${unit}`);
  if (result.c0 !== null) lines.push(`Initial concentration (C0) ${withUnit(result.c0, `${unit}/mL`)}`);
  if (result.cultureVolumeMl !== null) {
    lines.push(`Culture volume ${formatCount(result.cultureVolumeMl)} mL`);
  }
  if (result.blankDecayPercent !== null) lines.push(`Blank decay ${formatPercent(result.blankDecayPercent)}`);
  if (result.pePercent !== null) {
    lines.push(`Phytoremediation efficiency ${formatPercent(result.pePercent)}`);
  }
  if (result.balance) {
    lines.push(`Adjusted initial amount ${withUnit(result.balance.n0EffAdj, unit)}`);
    lines.push(`Water total ${withUnit(result.balance.nWater, unit)}`);
    lines.push(`Plant total ${withUnit(result.balance.nPlantTotal, unit)}`);
    lines.push(`Residual ${withUnit(result.balance.nResidual, unit)}`);
    lines.push(`Closure ${formatPercent(result.balance.closureFraction * 100)}`);
    lines.push(`Residual fraction ${formatPercent(result.balance.residualFraction * 100)}`);
  }
  if (result.tf !== null) lines.push(`Translocation factor ${formatCount(result.tf)}`);
  if (result.bcf !== null) lines.push(`Bioconcentration factor ${formatCount(result.bcf)}`);
  for (const part of result.parts.filter((row) => row.used)) {
    lines.push(
      `${part.name}: infective total ${formatCount(part.infective)}, inactivated total ${formatCount(part.inactivated)} and highly degraded particles ${formatCount(part.hdpTotal)}`,
    );
  }
  if (input.notes) lines.push(input.notes);
  return lines.join("\n");
}

function recompute() {
  renumberParts();
  latestInput = readForm();
  latestResult = computeSheet(latestInput);
  renderResults(latestInput, latestResult);
}

function loadJournal() {
  try {
    const raw = localStorage.getItem(STORAGE_KEY);
    const parsed = raw ? JSON.parse(raw) : [];
    return Array.isArray(parsed) ? parsed : [];
  } catch {
    return null;
  }
}

function writeJournal(records) {
  localStorage.setItem(STORAGE_KEY, JSON.stringify(records));
}

function refreshCounts() {
  const records = loadJournal();
  const text =
    records === null
      ? "Local storage is blocked in this browser."
      : records.length === 0
        ? "No records are stored in this browser."
        : records.length === 1
          ? "1 record is stored in this browser."
          : `${records.length} records are stored in this browser.`;
  $("browser-count").textContent = text;
  if ($("journal-count")) $("journal-count").textContent = text;
}

function renderJournal() {
  const records = loadJournal();
  const list = $("journal-list");
  const empty = $("journal-empty");
  list.replaceChildren();
  refreshCounts();
  if (records === null) {
    empty.hidden = false;
    empty.textContent = "Local storage is blocked in this browser. The calculation still runs.";
    return;
  }
  empty.hidden = records.length > 0;
  empty.textContent = "No records are stored in this browser. Save a record from the calculation sheet.";
  for (const record of records) {
    const item = document.createElement("li");
    const open = document.createElement("button");
    open.type = "button";
    open.className = "journal-open";
    open.textContent = record.title || "Untitled record";
    open.addEventListener("click", () => openRecord(record));
    const when = document.createElement("p");
    when.className = "journal-when";
    const parsed = new Date(record.savedAt);
    when.textContent = Number.isNaN(parsed.getTime())
      ? ""
      : parsed.toLocaleString(undefined, { dateStyle: "medium", timeStyle: "short" });
    const remove = document.createElement("button");
    remove.type = "button";
    remove.className = "secondary";
    remove.textContent = "Delete record";
    remove.addEventListener("click", () => deleteRecord(record.id, remove));
    item.append(open, when, remove);
    list.append(item);
  }
}

function openRecord(record) {
  currentId = record.id;
  fillForm(record.input);
  setStatus("Record opened.");
  location.hash = "#sheet";
}

function deleteRecord(id, button) {
  if (button.dataset.armed !== "yes") {
    button.dataset.armed = "yes";
    button.textContent = "Confirm deletion";
    window.setTimeout(() => {
      if (!button.isConnected) return;
      button.dataset.armed = "no";
      button.textContent = "Delete record";
    }, 4000);
    return;
  }
  const records = loadJournal();
  if (records === null) return;
  writeJournal(records.filter((record) => record.id !== id));
  if (currentId === id) currentId = null;
  renderJournal();
}

function saveSheet() {
  const records = loadJournal();
  if (records === null) {
    setStatus("Local storage is blocked. The calculation still runs. Export the record to keep a copy.");
    return;
  }
  const input = readForm();
  const id = currentId || crypto.randomUUID();
  currentId = id;
  const record = {
    id,
    savedAt: new Date().toISOString(),
    title: sheetTitle(input),
    input,
  };
  writeJournal([record, ...records.filter((item) => item.id !== id)]);
  refreshCounts();
  setStatus(`Record saved in this browser. ${sheetTitle(input)}.`);
}

function plainResult(result) {
  if (!result) return null;
  const balance = result.balance
    ? {
        n0EffAdj: result.balance.n0EffAdj,
        nWater: result.balance.nWater,
        nPlantTotal: result.balance.nPlantTotal,
        nResidual: result.balance.nResidual,
        closureFraction: result.balance.closureFraction,
        residualFraction: result.balance.residualFraction,
      }
    : null;
  return { ...result, balance };
}

function csvCell(value) {
  const text = value == null ? "" : String(value);
  if (/[",\n]/.test(text)) return `"${text.replaceAll('"', '""')}"`;
  return text;
}

function toCsv(input, result) {
  const rows = [["section", "key", "value"]];
  const add = (section, key, value) => rows.push([section, key, value ?? ""]);
  for (const [key, value] of Object.entries(input)) {
    if (key === "parts") continue;
    add("sheet", key, value);
  }
  input.parts.forEach((part, index) => {
    for (const [key, value] of Object.entries(part)) add(`part-${index + 1}`, key, value);
  });
  const plain = plainResult(result);
  if (plain) {
    for (const key of ["c0", "cultureVolumeMl", "blankDecayPercent", "pePercent", "tf", "bcf", "combinedPi", "mixturePi"]) {
      add("result", key, plain[key]);
    }
    if (plain.balance) {
      for (const [key, value] of Object.entries(plain.balance)) add("result", key, value);
    }
  }
  return rows.map((row) => row.map(csvCell).join(",")).join("\n");
}

function fileStamp() {
  return ($("date").value || today()).replaceAll("-", "");
}

function download(filename, text, type) {
  const blob = new Blob([text], { type });
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = filename;
  document.body.append(link);
  link.click();
  link.remove();
  window.setTimeout(() => URL.revokeObjectURL(url), 1500);
}

function exportJson() {
  recompute();
  download(
    `mass-balance-${fileStamp()}.json`,
    JSON.stringify({ input: latestInput, result: plainResult(latestResult) }, null, 2),
    "application/json",
  );
}

function exportCsv() {
  recompute();
  download(`mass-balance-${fileStamp()}.csv`, toCsv(latestInput, latestResult), "text/csv");
}

async function copyNotebook() {
  const text = $("notebook").textContent;
  try {
    await navigator.clipboard.writeText(text);
    setStatus("Results copied.");
  } catch {
    setStatus("Select the result text and copy it.");
  }
}

function loadExample(example) {
  currentId = null;
  clearArmed = false;
  $("clear-sheet").textContent = "Clear form";
  fillForm({ ...emptyInput(), ...example.input, date: today() });
  setStatus(example.note);
}

function clearSheet() {
  const button = $("clear-sheet");
  if (!clearArmed) {
    clearArmed = true;
    button.textContent = "Confirm clear form";
    window.setTimeout(() => {
      clearArmed = false;
      button.textContent = "Clear form";
    }, 4000);
    return;
  }
  clearArmed = false;
  button.textContent = "Clear form";
  currentId = null;
  fillForm(emptyInput());
  setStatus("Form cleared.");
}

function showView(name) {
  for (const view of ["sheet", "journal", "limits"]) {
    $(`view-${view}`).hidden = view !== name;
  }
  document.querySelectorAll("nav a").forEach((link) => {
    link.setAttribute("aria-current", link.dataset.view === name ? "page" : "false");
  });
  const bar = $("live-bar");
  if (name === "sheet") bar.removeAttribute("hidden");
  else bar.setAttribute("hidden", "");
  if (name === "journal") renderJournal();
}

function route(moveFocus) {
  const name = location.hash.replace("#", "");
  const view = name === "journal" || name === "limits" ? name : "sheet";
  showView(view);
  if (moveFocus) document.querySelector(`#view-${view} h1`)?.focus();
}

function renderNet() {
  $("net-status").textContent = navigator.onLine
    ? "Equations are calculated in this browser. A network connection is not required after the first visit."
    : "The browser is offline. The calculation sheet and stored records remain available.";
}

function bind() {
  const form = $("sheet-form");
  form.addEventListener("submit", (event) => event.preventDefault());
  form.addEventListener("input", (event) => {
    const part = event.target.closest?.("[data-part]");
    if (part && event.target.matches?.('[data-field="name"]')) syncLegend(part);
    recompute();
  });
  form.addEventListener("change", (event) => {
    if (event.target.id === "virus") applyVirusPreset();
    if (event.target.name === "blankMode") applyBlankMode();
    recompute();
  });
  $("add-part").addEventListener("click", () => {
    $("parts").append(createPart({ name: "" }));
    renumberParts();
    recompute();
  });
  $("parts").addEventListener("click", (event) => {
    const button = event.target.closest?.("[data-remove]");
    if (!button || button.disabled) return;
    button.closest("[data-part]")?.remove();
    renumberParts();
    recompute();
  });
  $("load-ms2").addEventListener("click", () => loadExample(ms2WorkedExample()));
  $("load-t4").addEventListener("click", () => loadExample(t4ReportedWater()));
  $("clear-sheet").addEventListener("click", clearSheet);
  $("save-sheet").addEventListener("click", saveSheet);
  $("export-json").addEventListener("click", exportJson);
  $("export-csv").addEventListener("click", exportCsv);
  $("copy-notebook").addEventListener("click", () => {
    void copyNotebook();
  });
  $("print-sheet").addEventListener("click", () => window.print());
  window.addEventListener("hashchange", () => route(true));
  window.addEventListener("online", renderNet);
  window.addEventListener("offline", renderNet);
}

bind();
fillForm(emptyInput());
refreshCounts();
renderNet();
route(false);

if ("serviceWorker" in navigator) {
  navigator.serviceWorker.register("./sw.js").catch(() => {});
}
