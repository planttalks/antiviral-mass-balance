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
  return title || "Untitled sheet";
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
  part.querySelector(".legend-name").textContent = name || "Part";
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
    item.textContent = error.message;
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

function resultLead(result) {
  if (result.errors.length) return "Fix the marked fields. Partial results stay on the sheet.";
  if (!result.balance) {
    const names = { spike: "the spike", blank: "the blank", water: "the water titer" };
    const missing = result.missing.map((key) => names[key]).filter(Boolean);
    if (missing.length) return `Enter ${joinWords(missing)} to close the balance.`;
    return "A harvested part is missing an infective total.";
  }
  const residual = result.balance.residualFraction;
  let lead = "The water and the parts match the adjusted spike.";
  if (residual > 0) lead = "Some of the adjusted spike is not in the water or the parts.";
  if (residual < 0) lead = "The water and the parts add up to more than the adjusted spike.";
  if (result.hdpIncluded) {
    return `${lead} HDp totals were added. Keep them in this closure only when they use the spike unit.`;
  }
  const showedHdp = result.parts.some((part) => part.hdpPerMl !== null);
  if (showedHdp) return `${lead} HDp is calculated and left out of the closure.`;
  return lead;
}

function renderResults(input, result) {
  const unit = input.unitSystem || "PFU";
  markErrors(result);
  $("result-lead").textContent = resultLead(result);
  $("spike-hint").textContent =
    result.c0 === null
      ? "C0 appears here after the stock and the volumes."
      : `C0 ${withUnit(result.c0, unit)}. Culture volume ${formatCount(result.cultureVolumeMl)} mL.`;
  $("blank-hint").textContent =
    result.blankDecayPercent === null
      ? "Blank decay is removed before the closure."
      : `Blank decay ${formatPercent(result.blankDecayPercent)}.`;
  $("water-hint").textContent =
    result.pePercent === null
      ? "Removal from the bulk water appears here."
      : `Removal from the water ${formatPercent(result.pePercent)}.`;

  const numbers = $("result-numbers");
  numbers.replaceChildren();
  addLine(numbers, "C0", withUnit(result.c0, `${unit}/mL`), "panel-c0");
  addLine(numbers, "Blank", formatPercent(result.blankDecayPercent), "panel-blank");
  addLine(numbers, "PE", formatPercent(result.pePercent), "panel-pe");
  if (result.balance) {
    addLine(numbers, "Adjusted spike", withUnit(result.balance.n0EffAdj, unit), "panel-n0");
    addLine(numbers, "Water", withUnit(result.balance.nWater, unit), "panel-water");
    addLine(numbers, "Plant", withUnit(result.balance.nPlantTotal, unit), "panel-plant");
    addLine(numbers, "Residual", withUnit(result.balance.nResidual, unit), "panel-residual-count");
    addLine(numbers, "Closure", formatPercent(result.balance.closureFraction * 100), "panel-closure");
    addLine(numbers, "Residual fraction", formatPercent(result.balance.residualFraction * 100), "panel-residual");
  } else {
    addLine(numbers, "Closure", "n/a", "panel-closure");
    addLine(numbers, "Residual fraction", "n/a", "panel-residual");
  }
  if (result.tf !== null) addLine(numbers, "TF", formatCount(result.tf), "panel-tf");
  if (result.bcf !== null) addLine(numbers, "BCF", formatCount(result.bcf), "panel-bcf");
  const weighted = result.parts.filter((part) => part.used && part.pi !== null && part.massG > 0);
  if (weighted.length >= 2 && result.combinedPi !== null) {
    addLine(numbers, "Eq. 4", withUnit(result.combinedPi, `${unit}/g`), "panel-eq4");
  }
  if (result.mixturePi !== null) {
    addLine(numbers, "Eq. 5", withUnit(result.mixturePi, `${unit}/g`), "panel-eq5");
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
    if (part.a !== null) lines.push(`A ${withUnit(part.a, `${unit}/mL`)}`);
    if (part.pi !== null) lines.push(`PI ${withUnit(part.pi, `${unit}/g`)}`);
    if (part.infective !== null) lines.push(`Infective ${withUnit(part.infective, unit)}`);
    if (part.inactivated !== null) lines.push(`Inactivated ${withUnit(part.inactivated, unit)}`);
    if (part.hdpPerMl !== null) lines.push(`HDp ${withUnit(part.hdpPerMl, "Vp/mL")}`);
    if (part.hdpTotal !== null) {
      lines.push(
        part.includedHdp
          ? `HDp total ${withUnit(part.hdpTotal, "Vp")} in the closure`
          : `HDp total ${withUnit(part.hdpTotal, "Vp")} left out`,
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
  lines.push(`Unit ${unit}`);
  if (result.c0 !== null) lines.push(`C0 ${withUnit(result.c0, `${unit}/mL`)}`);
  if (result.cultureVolumeMl !== null) lines.push(`Volume ${formatCount(result.cultureVolumeMl)} mL`);
  if (result.blankDecayPercent !== null) lines.push(`Blank ${formatPercent(result.blankDecayPercent)}`);
  if (result.pePercent !== null) lines.push(`PE ${formatPercent(result.pePercent)}`);
  if (result.balance) {
    lines.push(`Adjusted spike ${withUnit(result.balance.n0EffAdj, unit)}`);
    lines.push(`Water ${withUnit(result.balance.nWater, unit)}`);
    lines.push(`Plant ${withUnit(result.balance.nPlantTotal, unit)}`);
    lines.push(`Residual ${withUnit(result.balance.nResidual, unit)}`);
    lines.push(`Closure ${formatPercent(result.balance.closureFraction * 100)}`);
    lines.push(`Residual fraction ${formatPercent(result.balance.residualFraction * 100)}`);
  }
  if (result.tf !== null) lines.push(`TF ${formatCount(result.tf)}`);
  if (result.bcf !== null) lines.push(`BCF ${formatCount(result.bcf)}`);
  for (const part of result.parts.filter((row) => row.used)) {
    lines.push(
      `${part.name}: infective ${formatCount(part.infective)} inactivated ${formatCount(part.inactivated)} HDp ${formatCount(part.hdpTotal)}`,
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
      ? "This browser blocked local storage."
      : records.length === 0
        ? "No sheets on this browser."
        : `${records.length} ${records.length === 1 ? "sheet" : "sheets"} on this browser.`;
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
    empty.textContent = "This browser blocked local storage. The sheet still calculates.";
    return;
  }
  empty.hidden = records.length > 0;
  empty.textContent = "No sheets on this browser. Save one from the sheet.";
  for (const record of records) {
    const item = document.createElement("li");
    const open = document.createElement("button");
    open.type = "button";
    open.className = "journal-open";
    open.textContent = record.title || "Untitled sheet";
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
    remove.textContent = "Delete";
    remove.addEventListener("click", () => deleteRecord(record.id, remove));
    item.append(open, when, remove);
    list.append(item);
  }
}

function openRecord(record) {
  currentId = record.id;
  fillForm(record.input);
  setStatus("Opened a saved sheet.");
  location.hash = "#sheet";
}

function deleteRecord(id, button) {
  if (button.dataset.armed !== "yes") {
    button.dataset.armed = "yes";
    button.textContent = "Confirm delete";
    window.setTimeout(() => {
      if (!button.isConnected) return;
      button.dataset.armed = "no";
      button.textContent = "Delete";
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
    setStatus("This browser blocked local storage. Export the sheet instead.");
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
  setStatus(`Saved. ${sheetTitle(input)} is on this browser.`);
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
    setStatus("Copied the notebook lines.");
  } catch {
    setStatus("Select the notebook lines and copy them.");
  }
}

function loadExample(example) {
  currentId = null;
  clearArmed = false;
  $("clear-sheet").textContent = "Clear";
  fillForm({ ...emptyInput(), ...example.input, date: today() });
  setStatus(example.note);
}

function clearSheet() {
  const button = $("clear-sheet");
  if (!clearArmed) {
    clearArmed = true;
    button.textContent = "Confirm clear";
    window.setTimeout(() => {
      clearArmed = false;
      button.textContent = "Clear";
    }, 4000);
    return;
  }
  clearArmed = false;
  button.textContent = "Clear";
  currentId = null;
  fillForm(emptyInput());
  setStatus("Cleared.");
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
    ? "Equations run on this browser. A network is not required after the first load."
    : "Offline. The sheet and the journal still open.";
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
