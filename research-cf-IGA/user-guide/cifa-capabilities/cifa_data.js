// Shared constants, schema validation, data loading, filtering, and pagination helpers.

var DEFAULT_CARDS_PER_PAGE = 25;
var MIN_CARDS_PER_PAGE = 25;
var CARDS_PER_PAGE_STEP = 25;
var NAV_DEPTH_DISPLAY_LIMIT = 10;
var COREFORM_HOME_URL = "https://coreform.com/";
var COREFORM_LOGO_URL = "coreform-logo-wordmark-darkmode.svg";
var TITLE_LOGO_URL = "coreform-iga-icon.svg";
var CARDS_PER_PAGE_STORAGE_KEY = "coreform.cifaCapabilities.cardsPerPage.v1";
var ORDER_BY_CHANGES_FIRST_STORAGE_KEY = "coreform.cifaCapabilities.orderByChangesFirst.v1";
// Mirrors local_editor/schema.py CARD_ID_PATTERN.
var CARD_ID_PATTERN = /^[a-z0-9]+(?:-[a-z0-9]+)*$/;
var MONTH_YEAR_PATTERN = /^(January|February|March|April|May|June|July|August|September|October|November|December) [0-9]{4}$/;
var globalDismissHandlersBound = false;
var PLAN_MONTHS = [
  "January",
  "February",
  "March",
  "April",
  "May",
  "June",
  "July",
  "August",
  "September",
  "October",
  "November",
  "December"
];
var PLAN_YEARS = (function () {
  var years = [];
  var year;

  for (year = 2026; year <= 2050; year += 1) {
    years.push(String(year));
  }

  return years;
}());
var RELATIONSHIP_RECIPROCALS = {
  children: "parents",
  parents: "children"
};
var FILTER_FIELD_KEYS = ["tool", "status", "minAbaqusVersion", "minCoreformIgaVersion"];

function fullPatternMatch(pattern, value) {
  var match = typeof value === "string" ? pattern.exec(value) : null;
  return Boolean(match && match[0] === value);
}

function cardIdIsValid(value) {
  return fullPatternMatch(CARD_ID_PATTERN, value);
}

function escapeHtml(value) {
  return String(value).replace(/[&<>"']/g, function (character) {
    return {
      "&": "&amp;",
      "<": "&lt;",
      ">": "&gt;",
      '"': "&quot;",
      "'": "&#039;"
    }[character];
  });
}

function isObject(value) {
  return Boolean(value) && typeof value === "object" && !Array.isArray(value);
}

function validateFiniteJsonValues(value) {
  if (typeof value === "number" && !Number.isFinite(value)) {
    throw new Error("Capability sources must contain finite JSON values; an overflowing number cannot be saved.");
  }
  if (Array.isArray(value)) {
    value.forEach(validateFiniteJsonValues);
  } else if (isObject(value)) {
    Object.keys(value).forEach(function (key) { validateFiniteJsonValues(value[key]); });
  }
}

function parsePlanDate(value) {
  var parts = String(value || "").split(" ");

  if (parts.length === 2 && PLAN_MONTHS.indexOf(parts[0]) !== -1 && PLAN_YEARS.indexOf(parts[1]) !== -1) {
    return {
      month: parts[0],
      year: parts[1]
    };
  }

  return {
    month: "",
    year: ""
  };
}

function fetchJson(path) {
  if (!window.fetch) {
    return Promise.reject(new Error("This browser does not support fetch()."));
  }

  return fetch(path).then(function (response) {
    if (!response.ok) {
      throw new Error("HTTP " + response.status + " while fetching " + path + ".");
    }

    return response.json();
  });
}

function createEmptyRegistry() {
  return {
    version: 1,
    fields: [],
    fieldByKey: Object.create(null),
    requiredFields: [],
    tableFields: [],
    detailFields: [],
    editorFields: [],
    relationshipFields: [],
    textFields: []
  };
}

function configVisible(config) {
  return Boolean(config) && config.visible === true;
}

function validateBoolean(value, name) {
  if (value !== undefined && typeof value !== "boolean") {
    throw new Error(name + " must be a boolean.");
  }
}

function validateString(value, name) {
  if (typeof value !== "string" || value.trim() === "") {
    throw new Error(name + " must be a non-empty string.");
  }
}

function validateObjectConfig(value, name) {
  if (value === undefined) {
    return {};
  }
  if (!isObject(value)) {
    throw new Error(name + " must be an object.");
  }
  return value;
}

function validateStringArray(value, name) {
  var seen = Object.create(null);

  if (!Array.isArray(value) || !value.length) {
    throw new Error(name + " must be a non-empty array.");
  }

  value.forEach(function (item) {
    validateString(item, name + " item");
    if (seen[item]) {
      throw new Error(name + " contains duplicate value '" + item + "'.");
    }
    seen[item] = true;
  });
}

function validateKeyedLabelArray(value, name) {
  var seen = Object.create(null);

  if (!Array.isArray(value) || !value.length) {
    throw new Error(name + " must be a non-empty array.");
  }

  value.forEach(function (item, index) {
    var itemName = name + " item #" + (index + 1);

    if (!isObject(item)) {
      throw new Error(itemName + " must be an object.");
    }
    validateString(item.key, itemName + " key");
    validateString(item.label, itemName + " label");
    if (seen[item.key]) {
      throw new Error(name + " contains duplicate key '" + item.key + "'.");
    }
    seen[item.key] = true;
  });
}

function validateRegistry(value) {
  validateFiniteJsonValues(value);
  var registry = createEmptyRegistry();
  var allowedTypes = {
    string: true,
    text: true,
    enum: true,
    "card-list": true,
    "string-list": true,
    "url-list": true,
    plan: true,
    "boolean-map": true,
    "grouped-string-lists": true,
    object: true
  };
  var allowedRoles = { title: true, badge: true, summary: true, metadata: true, relationship: true };
  var allowedWidgets = {
    text: true,
    textarea: true,
    select: true,
    "card-picker": true,
    plan: true,
    "boolean-map": true,
    "grouped-string-lists": true,
    "production-basis": true
  };
  var hasIdField = false;

  if (!isObject(value)) {
    throw new Error("Expected capability-card-fields.json to contain an object.");
  }
  if (value.version !== 1) {
    throw new Error("capability-card-fields.json version must be 1.");
  }
  if (!Array.isArray(value.fields) || !value.fields.length) {
    throw new Error("capability-card-fields.json fields must be a non-empty array.");
  }

  value.fields.forEach(function (field, index) {
    var fieldName = "field #" + (index + 1);
    var table;
    var detail;
    var editor;

    if (!isObject(field)) {
      throw new Error(fieldName + " must be an object.");
    }

    validateString(field.key, fieldName + " key");
    validateString(field.label, "field '" + field.key + "' label");
    validateString(field.type, "field '" + field.key + "' type");
    validateBoolean(field.required, "field '" + field.key + "' required");
    validateBoolean(field.internalOnly, "field '" + field.key + "' internalOnly");
    validateBoolean(field.generated, "field '" + field.key + "' generated");

    if (registry.fieldByKey[field.key]) {
      throw new Error("Duplicate field key '" + field.key + "' in capability-card-fields.json.");
    }
    if (!Object.prototype.hasOwnProperty.call(allowedTypes, field.type)) {
      throw new Error("field '" + field.key + "' uses an unsupported type.");
    }
    if (field.type === "enum" || field.type === "plan") {
      validateStringArray(field.options, "field '" + field.key + "' options");
    }
    if (field.type === "boolean-map") {
      validateKeyedLabelArray(field.items, "field '" + field.key + "' items");
    }
    if (field.type === "grouped-string-lists") {
      validateKeyedLabelArray(field.groups, "field '" + field.key + "' groups");
    }

    table = validateObjectConfig(field.table, "field '" + field.key + "' table");
    validateBoolean(table.visible, "field '" + field.key + "' table.visible");
    if (configVisible(table)) {
      validateString(table.width, "field '" + field.key + "' table.width");
    }

    detail = validateObjectConfig(field.detail, "field '" + field.key + "' detail");
    validateBoolean(detail.visible, "field '" + field.key + "' detail.visible");
    if (configVisible(detail) && detail.role !== undefined) {
      validateString(detail.role, "field '" + field.key + "' detail.role");
      if (!Object.prototype.hasOwnProperty.call(allowedRoles, detail.role)) {
        throw new Error("field '" + field.key + "' detail.role must be title, badge, summary, metadata, or relationship.");
      }
    }

    editor = validateObjectConfig(field.editor, "field '" + field.key + "' editor");
    validateBoolean(editor.visible, "field '" + field.key + "' editor.visible");
    validateBoolean(editor.readOnlyOnEdit, "field '" + field.key + "' editor.readOnlyOnEdit");
    if (configVisible(editor)) {
      validateString(editor.widget, "field '" + field.key + "' editor.widget");
    if (!Object.prototype.hasOwnProperty.call(allowedWidgets, editor.widget)) {
        throw new Error("field '" + field.key + "' editor.widget must be text, textarea, select, card-picker, plan, boolean-map, grouped-string-lists, or production-basis.");
    }
      if (editor.widget === "select" && field.type !== "enum") {
        throw new Error("field '" + field.key + "' uses a select editor but is not an enum.");
      }
      if (editor.widget === "card-picker" && field.type !== "card-list") {
        throw new Error("field '" + field.key + "' uses a card-picker editor but is not a card-list.");
      }
      if (editor.widget === "plan" && field.type !== "plan") {
        throw new Error("field '" + field.key + "' uses a plan editor but is not a plan field.");
      }
      if (editor.widget === "boolean-map" && field.type !== "boolean-map") {
        throw new Error("field '" + field.key + "' uses a boolean-map editor but is not a boolean-map field.");
      }
      if (editor.widget === "grouped-string-lists" && field.type !== "grouped-string-lists") {
        throw new Error("field '" + field.key + "' uses a grouped-string-lists editor but is not a grouped-string-lists field.");
      }
      if (editor.widget === "production-basis" && (field.type !== "object" || field.key !== "productionBasis")) {
        throw new Error("The production-basis editor requires the productionBasis object field.");
      }
    }

    registry.fieldByKey[field.key] = field;
    registry.fields.push(field);
    if (field.required === true) {
      registry.requiredFields.push(field);
    }
    if (configVisible(table)) {
      registry.tableFields.push(field);
    }
    if (configVisible(detail)) {
      registry.detailFields.push(field);
    }
    if (configVisible(editor)) {
      registry.editorFields.push(field);
    }
    if (field.type === "card-list") {
      registry.relationshipFields.push(field);
    }
    if (field.type === "text" || field.type === "string") {
      registry.textFields.push(field);
    }
    if (field.key === "id") {
      hasIdField = true;
    }
  });

  if (!hasIdField) {
    throw new Error("capability-card-fields.json must define an 'id' field.");
  }

  return registry;
}

function parsedInlineLinkTarget(source) {
  var trimmed = String(source).trim();
  var labelMatch = trimmed.match(/^(.+)\s+<([^>]+)>$/);
  var label = "";
  var target = trimmed;

  if (labelMatch) {
    label = labelMatch[1].trim();
    target = labelMatch[2].trim();
    if (label === "" || target === "") {
      return null;
    }
  }

  if (validUrl(target)) {
    return {
      type: "external",
      url: target,
      label: label || target
    };
  }

  if (cardIdIsValid(target)) {
    return {
      type: "card",
      id: target,
      label: label
    };
  }

  return null;
}

function textCardLinkIds(value, fieldLabel) {
  var matches = [];

  String(value).replace(/\[\[([^\]]+)\]\]/g, function (match, source) {
    var target = parsedInlineLinkTarget(source);

    if (!target) {
      throw new Error("'" + fieldLabel + "' has invalid link syntax '" + match + "'.");
    }
    if (target.type === "card") {
      matches.push(target.id);
    }
    return match;
  });

  return matches;
}

function validateCardListValue(value, field, card, cardName) {
  var seen = Object.create(null);

  if (!Array.isArray(value)) {
    throw new Error(cardName + " field '" + field.key + "' must be an array of relationship entries.");
  }

  value.forEach(function (item) {
    var entry = relationshipEntryFromValue(item);
    var cardId = entry.id;

    if (typeof item !== "string" && !isObject(item)) {
      throw new Error(cardName + " field '" + field.key + "' contains an invalid relationship entry.");
    }
    if (isObject(item) && typeof item.supported !== "boolean") {
      throw new Error(cardName + " field '" + field.key + "' entry for '" + cardId + "' must include boolean supported.");
    }
    if (isObject(item) && (Object.keys(item).length !== 2 || !Object.prototype.hasOwnProperty.call(item, "id"))) {
      throw new Error(cardName + " relationship entries require exactly id and supported.");
    }
    if (!cardIdIsValid(cardId)) {
      throw new Error(cardName + " field '" + field.key + "' contains an invalid card ID.");
    }
    if (cardId === card.id) {
      throw new Error(cardName + " field '" + field.key + "' cannot reference itself.");
    }
    if (seen[cardId]) {
      throw new Error(cardName + " field '" + field.key + "' contains duplicate card ID '" + cardId + "'.");
    }
    seen[cardId] = true;
  });
}

function validUrl(value) {
  var url;

  if (typeof value !== "string" || !/^https?:\/\//i.test(value) || /[\s\x00-\x1f\\]/.test(value)) {
    return false;
  }

  try {
    url = new URL(value);
  } catch (error) {
    return false;
  }

  return (url.protocol === "http:" || url.protocol === "https:") && Boolean(url.hostname) && !url.username && !url.password;
}

function validManualPage(value) {
  return validUrl(value) || (typeof value === "string" && !/[\s\x00-\x1f\\]/.test(value) &&
    /^\.\.\/verification-manual\/(?:[A-Za-z0-9_-]+\/)*[A-Za-z0-9_-]+\.html(?:#[A-Za-z0-9_.:-]+)?$/.test(value));
}

function validateStringListValue(value, field, cardName) {
  if (!Array.isArray(value)) {
    throw new Error(cardName + " field '" + field.key + "' must be an array of strings.");
  }
  value.forEach(function (item) {
    if (typeof item !== "string" || item.trim() === "") {
      throw new Error(cardName + " field '" + field.key + "' must contain only non-empty strings.");
    }
  });
}

function validateUrlListValue(value, field, cardName) {
  if (!Array.isArray(value)) {
    throw new Error(cardName + " field '" + field.key + "' must be an array of URLs.");
  }
  value.forEach(function (item) {
    if (!(field.key === "verificationTests" ? validManualPage(item) : validUrl(item))) {
      throw new Error(cardName + " field '" + field.key + "' contains an invalid URL.");
    }
  });
}

function validatePlanValue(value, field, cardName) {
  if (!isObject(value)) {
    throw new Error(cardName + " field '" + field.key + "' must be a plan object.");
  }
  requireKnownKeys(value, ["stage", "date"], cardName + " plan");
  if (field.options.indexOf(value.stage) === -1) {
    throw new Error(cardName + " field '" + field.key + "' has an invalid stage.");
  }
  if (!fullPatternMatch(MONTH_YEAR_PATTERN, value.date)) {
    throw new Error(cardName + " field '" + field.key + "' date must use 'Month YYYY'.");
  }
}

function validateBooleanMapValue(value, field, cardName) {
  var allowed = Object.create(null);

  if (!isObject(value)) {
    throw new Error(cardName + " field '" + field.key + "' must be an object of booleans.");
  }
  field.items.forEach(function (item) {
    allowed[item.key] = true;
  });
  Object.keys(value).forEach(function (key) {
    if (!allowed[key]) {
      throw new Error(cardName + " field '" + field.key + "' contains unknown key '" + key + "'.");
    }
    if (typeof value[key] !== "boolean") {
      throw new Error(cardName + " field '" + field.key + "' value for '" + key + "' must be boolean.");
    }
  });
}

function validateGroupedStringListsValue(value, field, cardName) {
  var allowed = Object.create(null);

  if (!isObject(value)) {
    throw new Error(cardName + " field '" + field.key + "' must be an object of string arrays.");
  }
  field.groups.forEach(function (group) {
    allowed[group.key] = true;
  });
  Object.keys(value).forEach(function (key) {
    if (!allowed[key]) {
      throw new Error(cardName + " field '" + field.key + "' contains unknown group '" + key + "'.");
    }
    if (!Array.isArray(value[key])) {
      throw new Error(cardName + " field '" + field.key + "' group '" + key + "' must be an array.");
    }
    value[key].forEach(function (item) {
      if (typeof item !== "string" || item.trim() === "") {
        throw new Error(cardName + " field '" + field.key + "' group '" + key + "' must contain only non-empty strings.");
      }
    });
  });
}

function fieldHasRequiredValue(value, field) {
  if (value === undefined) {
    return false;
  }
  if (field.type === "card-list" || field.type === "string-list" || field.type === "url-list") {
    return Array.isArray(value) && value.length > 0;
  }
  if (field.type === "plan" || field.type === "boolean-map" || field.type === "grouped-string-lists" || field.type === "object") {
    return isObject(value) && Object.keys(value).length > 0;
  }
  return typeof value === "string" && value.trim() !== "";
}

function validateCardFieldValue(value, field, card, cardName) {
  if (field.type === "card-list") {
    validateCardListValue(value, field, card, cardName);
    return;
  }
  if (field.type === "string-list") {
    validateStringListValue(value, field, cardName);
    return;
  }
  if (field.type === "url-list") {
    validateUrlListValue(value, field, cardName);
    return;
  }
  if (field.type === "plan") {
    validatePlanValue(value, field, cardName);
    return;
  }
  if (field.type === "boolean-map") {
    validateBooleanMapValue(value, field, cardName);
    return;
  }
  if (field.type === "grouped-string-lists") {
    validateGroupedStringListsValue(value, field, cardName);
    return;
  }
  if (field.type === "object") {
    if (!isObject(value)) {
      throw new Error(cardName + " field '" + field.key + "' must be an object.");
    }
    return;
  }
  if (typeof value !== "string") {
    throw new Error(cardName + " field '" + field.key + "' must be a string.");
  }
  if (field.type === "enum" && field.options.indexOf(value) === -1) {
    throw new Error(cardName + " field '" + field.key + "' value '" + value + "' is not in the registry options.");
  }
  if (field.type === "text" || field.type === "string") {
    textCardLinkIds(value, field.label);
  }
}

function validateRelationshipSides(card, cardName) {
  var childIds = Object.create(null);

  cardListValue(card, "children").forEach(function (entry) {
    childIds[entry.id] = true;
  });
  cardListValue(card, "parents").forEach(function (entry) {
    if (childIds[entry.id]) {
      throw new Error(cardName + " cannot list card ID '" + entry.id + "' as both a child and a parent.");
    }
  });
}

function reciprocalRelationshipKey(key) {
  return RELATIONSHIP_RECIPROCALS[key] || "";
}

function validateCardReferenceTargets(card, cardName, registry, ids, cardsById) {
  registry.relationshipFields.forEach(function (field) {
    cardListValue(card, field.key).forEach(function (entry) {
      var cardId = entry.id;
      var reciprocalKey;
      var targetCard;
      var reciprocalEntry;

      if (!ids[cardId]) {
        throw new Error(cardName + " field '" + field.key + "' references missing card ID '" + cardId + "'.");
      }
      reciprocalKey = reciprocalRelationshipKey(field.key);
      targetCard = cardsById && cardsById[cardId];
      reciprocalEntry = reciprocalKey && targetCard ? relationshipEntryFor(cardListValue(targetCard, reciprocalKey), card.id) : null;
      if (reciprocalKey && targetCard && !reciprocalEntry) {
        throw new Error(
          cardName + " field '" + field.key + "' references card ID '" + cardId +
          "', but that card does not reciprocate with field '" + reciprocalKey + "'."
        );
      }
      if (reciprocalEntry && reciprocalEntry.supported !== entry.supported) {
        throw new Error(
          cardName + " field '" + field.key + "' references card ID '" + cardId +
          "', but that card has a different supported value in field '" + reciprocalKey + "'."
        );
      }
    });
  });

  registry.textFields.forEach(function (field) {
    var value = card[field.key];

    if (typeof value === "string") {
      textCardLinkIds(value, field.label).forEach(function (cardId) {
        if (!ids[cardId]) {
          throw new Error(cardName + " field '" + field.key + "' links to missing card ID '" + cardId + "'.");
        }
      });
    }
  });
}

function validateCards(value, registry, options) {
  validateFiniteJsonValues(value);
  var ids = Object.create(null);
  var cardsById = Object.create(null);

  if (!Array.isArray(value)) {
    throw new Error("Expected capability-cards.json to contain an array.");
  }

  value.forEach(function (card, index) {
    var cardName = "card #" + (index + 1);

    if (!isObject(card)) {
      throw new Error(cardName + " must be a JSON object.");
    }

    Object.keys(card).forEach(function (key) {
      if (!registry.fieldByKey[key]) {
        throw new Error(cardName + " has unknown field '" + key + "'. Add it to capability-card-fields.json first.");
      }
      if (options && options.source === true && registry.fieldByKey[key].generated === true) {
        throw new Error("Source field '" + key + "' is generated; remove it and regenerate evidence.");
      }
    });

    registry.requiredFields.forEach(function (field) {
      var fieldValue = card[field.key];

      if (!fieldHasRequiredValue(fieldValue, field)) {
        throw new Error(cardName + " is missing required field '" + field.key + "'.");
      }
    });

    registry.fields.forEach(function (field) {
      var value = card[field.key];

      if (value === undefined) {
        return;
      }
      validateCardFieldValue(value, field, card, cardName);
    });

    if (!cardIdIsValid(card.id)) {
      throw new Error("Card ID must be a lowercase hyphenated stable ID, for example 'analysis-linear-static'.");
    }
    if (ids[card.id]) {
      throw new Error("Duplicate card id '" + card.id + "' in capability-cards.json.");
    }

    validateRelationshipSides(card, cardName);
    if (options && options.source === true) {
      validateProductionBasis(card);
    }
    ids[card.id] = true;
    cardsById[card.id] = card;
  });

  if (!options || options.checkGraph !== false) {
    value.forEach(function (card, index) {
      validateCardReferenceTargets(card, "card #" + (index + 1), registry, ids, cardsById);
    });
    validateAcyclicRelationships(cardsById);
  }

  return value;
}

function readApiJson(response) {
  return response.text().then(function (text) {
    var data = {};

    if (text) {
      try {
        data = JSON.parse(text);
      } catch (error) {
        throw new Error("The local editor API returned invalid JSON.");
      }
    }

    if (!response.ok || data.success === false) {
      throw new Error(data.message || ("Local editor API failed with HTTP " + response.status + "."));
    }

    return data;
  });
}

function loadInternalState() {
  if (!window.fetch) {
    return Promise.reject(new Error("This browser does not support fetch()."));
  }

  return fetch("/api/cifa-capabilities/state", {
    credentials: "same-origin",
    headers: {
      "Accept": "application/json"
    }
  }).then(readApiJson).then(function (data) {
    var registry = validateRegistry(data.registry);

    return {
      cards: validateCards(data.cards || [], registry, editorCollectionOptions(data)),
      registry: registry,
      changeMetadataById: isObject(data.changeMetadataById) ? data.changeMetadataById : {},
      gitStatus: isObject(data.gitStatus) ? data.gitStatus : null,
      evidenceStatus: isObject(data.evidenceStatus) ? data.evidenceStatus : null,
      eligibilityById: isObject(data.eligibilityById) ? data.eligibilityById : {},
      changedFiles: Array.isArray(data.changedFiles) ? data.changedFiles : [],
      collectionVersion: typeof data.collectionVersion === "string" ? data.collectionVersion : null,
      diagnostics: Array.isArray(data.diagnostics) ? data.diagnostics : [],
      writesAllowed: data.writesAllowed !== false
    };
  });
}

function editorCollectionOptions(data) {
  var repairing = data.writesAllowed === false && Array.isArray(data.diagnostics) && data.diagnostics.length > 0;
  return { source: false, checkGraph: !repairing };
}

function loadPublicData() {
  return Promise.all([
    fetchJson("capability-cards.json"),
    fetchJson("capability-card-fields.json")
  ]).then(function (payloads) {
    var registry = validateRegistry(payloads[1]);

    return {
      cards: validateCards(payloads[0], registry),
      registry: registry,
      changeMetadataById: {}
    };
  });
}

function loadData(mode) {
  return mode === "internal" ? loadInternalState() : loadPublicData();
}

function cardById(cards, id) {
  for (var index = 0; index < cards.length; index += 1) {
    if (cards[index].id === id) {
      return cards[index];
    }
  }
  return null;
}

function cardHref(cardId) {
  return "?card=" + encodeURIComponent(cardId);
}

function updateCardUrl(cardId) {
  var url;

  if (!window.history || !window.history.replaceState) {
    return;
  }

  try {
    url = new URL(window.location.href);
    if (cardId) {
      url.searchParams.set("card", cardId);
    } else {
      url.searchParams.delete("card");
    }
    window.history.replaceState({}, "", url.pathname + url.search + url.hash);
  } catch (error) {
    return;
  }
}

function cardIdFromLocation() {
  var match;

  try {
    if (window.URLSearchParams) {
      return new URLSearchParams(window.location.search).get("card") || "";
    }

    match = window.location.search.match(/[?&]card=([^&]+)/);
    return match ? decodeURIComponent(match[1].replace(/\+/g, " ")) : "";
  } catch (error) {
    return "";
  }
}

function plainLeftClick(event) {
  return event.button === 0 && !event.metaKey && !event.ctrlKey && !event.shiftKey && !event.altKey;
}

function fieldValue(card, field) {
  if (!card || card[field.key] === undefined) {
    return "";
  }
  return card[field.key];
}

function displayFieldValue(card, field) {
  var value = fieldValue(card, field);

  if (field.type === "card-list") {
    return cardListValue(card, field.key).map(function (entry) {
      return relationshipDisplayLabel(entry);
    }).join(", ");
  }
  if (field.type === "string-list" || field.type === "url-list") {
    return Array.isArray(value) ? value.join(", ") : "";
  }
  if (field.type === "plan") {
    return isObject(value) ? [value.stage, value.date].filter(Boolean).join(" ") : "";
  }
  if (field.type === "boolean-map") {
    return isObject(value) ? field.items.filter(function (item) {
      return value[item.key] === true;
    }).map(function (item) {
      return item.label;
    }).join(", ") : "";
  }
  if (field.type === "grouped-string-lists") {
    return isObject(value) ? field.groups.map(function (group) {
      return Array.isArray(value[group.key]) && value[group.key].length ? group.label + ": " + value[group.key].join(", ") : "";
    }).filter(Boolean).join("; ") : "";
  }

  if (Array.isArray(value)) {
    return value.join(", ");
  }

  return value;
}

function relationshipEntryFromValue(value) {
  if (typeof value === "string") {
    return {
      id: value,
      supported: true
    };
  }
  if (isObject(value)) {
    return {
      id: value.id,
      supported: value.supported
    };
  }
  return {
    id: "",
    supported: true
  };
}

function normalizeRelationshipEntries(value) {
  if (!Array.isArray(value)) {
    return [];
  }

  return value.map(relationshipEntryFromValue);
}

function cardListValue(card, key) {
  var value = card && card[key];

  return normalizeRelationshipEntries(value);
}

function relationshipDisplayLabel(entry) {
  return entry.id + " (" + (entry.supported ? "supported" : "not supported") + ")";
}

function relationshipEntryFor(entries, cardId) {
  for (var index = 0; index < entries.length; index += 1) {
    if (entries[index].id === cardId) {
      return entries[index];
    }
  }
  return null;
}

function cloneCard(card) {
  var copy = {};

  Object.keys(card || {}).forEach(function (key) {
    copy[key] = Array.isArray(card[key]) ? card[key].map(function (item) {
      return isObject(item) ? Object.assign({}, item) : item;
    }) : card[key];
  });

  return copy;
}

function cloneCards(cards) {
  return cards.map(function (card) {
    return cloneCard(card);
  });
}

function loadStoredCardsPerPagePreference() {
  var storedText;

  try {
    storedText = window.localStorage.getItem(CARDS_PER_PAGE_STORAGE_KEY);
  } catch (error) {
    return DEFAULT_CARDS_PER_PAGE;
  }

  return storedText ? normalizeCardsPerPage(storedText) : DEFAULT_CARDS_PER_PAGE;
}

function saveStoredCardsPerPagePreference(cardsPerPage) {
  try {
    window.localStorage.setItem(CARDS_PER_PAGE_STORAGE_KEY, String(normalizeCardsPerPage(cardsPerPage)));
  } catch (error) {
    return false;
  }

  return true;
}

function loadStoredOrderByChangesFirstPreference() {
  var storedText;

  try {
    storedText = window.localStorage.getItem(ORDER_BY_CHANGES_FIRST_STORAGE_KEY);
  } catch (error) {
    return false;
  }

  return storedText === "true";
}

function saveStoredOrderByChangesFirstPreference(enabled) {
  try {
    window.localStorage.setItem(ORDER_BY_CHANGES_FIRST_STORAGE_KEY, enabled ? "true" : "false");
  } catch (error) {
    return false;
  }

  return true;
}

function cardIdsForCards(cards) {
  var ids = Object.create(null);

  cards.forEach(function (card) {
    ids[card.id] = true;
  });

  return ids;
}

function cardMapForCards(cards) {
  var map = Object.create(null);

  cards.forEach(function (card) {
    map[card.id] = card;
  });

  return map;
}

function requireKnownKeys(value, allowed, label) {
  Object.keys(value).forEach(function (key) {
    if (allowed.indexOf(key) === -1) {
      throw new Error(label + " contains unknown key '" + key + "'.");
    }
  });
}

function validateArtifactGlob(value, label) {
  if (typeof value !== "string" || !value.trim() || value[0] === "/" || /[\\\x00:]/.test(value) || value.split("/").indexOf("..") !== -1) {
    throw new Error(label + " must be a workspace-relative glob without parent traversal.");
  }
}

// Python owns regex syntax; Create/Save validates patterns on the server before any write.
function validateEvidenceContract(contract) {
  var matcherKeys = {
    artifactText: ["glob", "contains", "keywords"], journal: ["glob", "contains", "keywords"],
    sourceInp: ["glob", "contains", "keywords"], translatedInp: ["glob", "contains", "keywords"],
    performanceSummary: ["glob", "matches"], checks: ["patterns", "consumes"],
    keywordReport: ["glob", "keywords", "supportCaseIds"], sqlite: ["glob", "tables"]
  };
  var caseKeys = ["name", "ctestPattern", "slug", "artifacts", "resultFields", "manualPage"].concat(Object.keys(matcherKeys));
  var names = Object.create(null);
  if (!isObject(contract)) {
    throw new Error("evidenceContract must be an object.");
  }
  requireKnownKeys(contract, ["minimumLevel", "anyOf"], "evidenceContract");
  if (contract.minimumLevel !== undefined && ["Exercised", "Numerically validated"].indexOf(contract.minimumLevel) === -1) {
    throw new Error("evidenceContract.minimumLevel must be Exercised or Numerically validated.");
  }
  if (!Array.isArray(contract.anyOf) || !contract.anyOf.length) {
    throw new Error("evidenceContract.anyOf must be a nonempty array.");
  }
  contract.anyOf.forEach(function (item) {
    if (!isObject(item)) {
      throw new Error("Each evidence case must be an object.");
    }
    requireKnownKeys(item, caseKeys, "Evidence case");
    validateString(item.name, "Evidence case name");
    if (names[item.name]) {
      throw new Error("Evidence case names must be unique.");
    }
    names[item.name] = true;
    validateString(item.ctestPattern, "Evidence case ctestPattern");
    if (item.manualPage !== undefined && !validManualPage(item.manualPage)) {
      throw new Error("Evidence case manualPage must be an HTTP(S) URL or ../verification-manual/...html link.");
    }
    if (item.slug !== undefined && !isObject(item.slug)) {
      throw new Error("Evidence case slug must be an object.");
    }
    ["artifacts", "resultFields"].forEach(function (key) {
      if (item[key] !== undefined) {
        validateStringListValue(item[key], { key: key }, "Evidence case");
      }
    });
    (item.artifacts || []).forEach(function (glob) { validateArtifactGlob(glob, "Artifact"); });
    Object.keys(matcherKeys).forEach(function (key) {
      if (item[key] === undefined) {
        return;
      }
      var spec = item[key];
      if (!isObject(spec)) {
        throw new Error("Evidence case " + key + " must be an object.");
      }
      requireKnownKeys(spec, matcherKeys[key], "Evidence case " + key);
      Object.keys(spec).forEach(function (part) {
        if (part === "glob") {
          validateArtifactGlob(spec[part], key + ".glob");
        } else if (part === "matches") {
          if (!isObject(spec[part])) {
            throw new Error(key + ".matches must be an object.");
          }
        } else {
          validateStringListValue(spec[part], { key: part }, key);
        }
      });
    });
  });
}

function validateProductionBasis(card) {
  var basis = card.productionBasis;
  if (Object.prototype.hasOwnProperty.call(card, "manualEvidence")) {
    throw new Error("Replace legacy manualEvidence with an explicit justified Production override.");
  }
  if (basis === undefined) {
    if (card.evidenceContract !== undefined) {
      throw new Error("evidenceContract requires productionBasis.kind=verification.");
    }
    if (card.status === "Production" && !cardListValue(card, "children").some(function (entry) { return entry.supported === true; })) {
      throw new Error("Production requires a verification or override basis, or supported children.");
    }
    return;
  }
  if (!isObject(basis) || ["verification", "override"].indexOf(basis.kind) === -1) {
    throw new Error("productionBasis.kind must be verification or override.");
  }
  if (basis.kind === "verification") {
    requireKnownKeys(basis, ["kind"], "Verification basis");
    validateEvidenceContract(card.evidenceContract);
  } else {
    requireKnownKeys(basis, ["kind", "justification", "references"], "Override basis");
    if (card.evidenceContract !== undefined) {
      throw new Error("Override cannot also declare evidenceContract.");
    }
    validateString(basis.justification, "Override justification");
    if (basis.references !== undefined) {
      validateUrlListValue(basis.references, { key: "references" }, "Override");
    }
  }
}

function validateAcyclicRelationships(cardsById) {
  var visited = Object.create(null);
  var active = Object.create(null);
  Object.keys(cardsById).forEach(function (start) {
    var stack = [[start, false]];
    var trail = [];
    while (stack.length) {
      var next = stack.pop();
      var id = next[0];
      if (next[1]) {
        delete active[id];
        visited[id] = true;
        trail.pop();
      } else if (active[id]) {
        throw new Error("Capability relationship cycle: " + trail.slice(trail.indexOf(id)).concat([id]).join(" -> "));
      } else if (!visited[id]) {
        active[id] = true;
        trail.push(id);
        stack.push([id, true]);
        cardListValue(cardsById[id], "children").slice().reverse().forEach(function (entry) {
          stack.push([entry.id, false]);
        });
      }
    }
  });
}

function validateCardForEditor(card, registry, cards, options) {
  validateFiniteJsonValues(card);
  var ids = cards ? cardIdsForCards(cards) : null;
  var settings = options || {};
  var cardsById = cards && settings.checkReciprocal !== false ? cardMapForCards(cards) : null;

  if (!isObject(card)) {
    throw new Error("Card data must be a JSON object.");
  }

  Object.keys(card).forEach(function (key) {
    if (!registry.fieldByKey[key]) {
      throw new Error("Unknown field '" + key + "'. Add it to capability-card-fields.json first.");
    }
    if (settings.source !== false && registry.fieldByKey[key].generated === true) {
      throw new Error("Source field '" + key + "' is generated; remove it and regenerate evidence.");
    }
  });

  registry.requiredFields.forEach(function (field) {
    var value = card[field.key];

    if (!fieldHasRequiredValue(value, field)) {
      throw new Error("'" + field.label + "' is required.");
    }
  });

  registry.fields.forEach(function (field) {
    var value = card[field.key];

    if (value === undefined) {
      return;
    }
    validateCardFieldValue(value, field, card, "Card '" + (card.id || "(new)") + "'");
  });

  if (!cardIdIsValid(card.id)) {
    throw new Error("Card ID must be a lowercase hyphenated stable ID: lowercase letters and digits separated by single hyphens, for example 'analysis-linear-static'.");
  }

  validateRelationshipSides(card, "Card '" + card.id + "'");
  if (settings.source !== false) {
    validateProductionBasis(card);
  }

  if (ids) {
    validateCardReferenceTargets(card, "Card '" + card.id + "'", registry, ids, cardsById);
  }
}

function normalizeCardsPerPage(value) {
  var numericValue = parseInt(value, 10);
  var steppedValue;

  if (!numericValue) {
    return DEFAULT_CARDS_PER_PAGE;
  }

  steppedValue = Math.round(numericValue / CARDS_PER_PAGE_STEP) * CARDS_PER_PAGE_STEP;

  if (steppedValue < MIN_CARDS_PER_PAGE) {
    return MIN_CARDS_PER_PAGE;
  }

  return steppedValue;
}

function pageForCard(cards, id, cardsPerPage) {
  var index;
  var pageSize = normalizeCardsPerPage(cardsPerPage);

  for (index = 0; index < cards.length; index += 1) {
    if (cards[index].id === id) {
      return Math.floor(index / pageSize) + 1;
    }
  }

  return 1;
}

function sanitizeCssTrack(value) {
  var track = typeof value === "string" ? value.trim() : "";

  if (!track || !/^[A-Za-z0-9 .,%()+-]+$/.test(track)) {
    return "1fr";
  }

  return track;
}

function tableGridTemplate(registry) {
  return registry.tableFields.map(function (field) {
    return "minmax(0, " + sanitizeCssTrack(field.table.width) + ")";
  }).join(" ");
}

function tableGridStyle(registry) {
  var template = tableGridTemplate(registry);

  return template ? ' style="grid-template-columns: ' + escapeHtml(template) + '"' : "";
}

function fieldByDetailRole(registry, role) {
  for (var index = 0; index < registry.detailFields.length; index += 1) {
    if (registry.detailFields[index].detail.role === role) {
      return registry.detailFields[index];
    }
  }
  return null;
}

function changeMetadataForCard(changeMetadataById, cardId) {
  if (!changeMetadataById || !isObject(changeMetadataById)) {
    return null;
  }

  return changeMetadataById[cardId] || null;
}

function changeMetadataBadge(metadata) {
  if (!metadata) {
    return "";
  }

  if (metadata.type === "create") {
    return '<span class="cm-change-pill cm-change-pill-new">New</span>';
  }
  if (metadata.type === "edit") {
    return '<span class="cm-change-pill cm-change-pill-edited">Edited</span>';
  }

  return "";
}

function changeSortRank(metadata) {
  if (!metadata) {
    return 2;
  }
  if (metadata.type === "create") {
    return 0;
  }
  if (metadata.type === "edit") {
    return 1;
  }

  return 2;
}

function cardSortTitle(card, registry) {
  var titleField = fieldByDetailRole(registry, "title");
  var title = titleField ? fieldValue(card, titleField) : "";

  return String(title || card.id || "").toLowerCase();
}

function cardMatchesSearch(card, registry, query) {
  var normalizedQuery = String(query || "").trim().toLowerCase();

  if (!normalizedQuery) {
    return true;
  }

  return cardSortTitle(card, registry).indexOf(normalizedQuery) !== -1;
}

function createEmptyFilters() {
  return {
    tool: "",
    status: "",
    minAbaqusVersion: "",
    minCoreformIgaVersion: ""
  };
}

function cardMatchesFieldFilters(card, registry, filters) {
  var activeFilters = filters || createEmptyFilters();

  return FILTER_FIELD_KEYS.every(function (key) {
    var value = String(activeFilters[key] || "");
    var field = registry.fieldByKey[key];

    if (!value) {
      return true;
    }

    return field ? displayFieldValue(card, field) === value : true;
  });
}

function hasActiveCardFilters(state) {
  if (String(state.searchQuery || "").trim()) {
    return true;
  }

  return FILTER_FIELD_KEYS.some(function (key) {
    return Boolean(state.filters && state.filters[key]);
  });
}

function displayCardsForState(state, registry) {
  var cards = state.cards.filter(function (card) {
    return cardMatchesSearch(card, registry, state.searchQuery) && cardMatchesFieldFilters(card, registry, state.filters);
  });

  if (!state.orderByChangesFirst) {
    return cards;
  }

  return cards.slice().sort(function (left, right) {
    var leftRank = changeSortRank(changeMetadataForCard(state.changeMetadataById, left.id));
    var rightRank = changeSortRank(changeMetadataForCard(state.changeMetadataById, right.id));
    var leftTitle;
    var rightTitle;

    if (leftRank !== rightRank) {
      return leftRank - rightRank;
    }

    leftTitle = cardSortTitle(left, registry);
    rightTitle = cardSortTitle(right, registry);
    if (leftTitle < rightTitle) {
      return -1;
    }
    if (leftTitle > rightTitle) {
      return 1;
    }

    return String(left.id || "").localeCompare(String(right.id || ""));
  });
}
