// Local-editor API calls, state transitions, event binding, and application startup.

function sendEditorRequest(method, path, payload) {
  var options = {
    method: method,
    credentials: "same-origin",
    headers: {
      "Accept": "application/json"
    }
  };

  if (payload !== undefined) {
    options.headers["Content-Type"] = "application/json";
    options.body = JSON.stringify(payload);
  }

  return fetch(path, options).then(readApiJson);
}

function applyEditorState(root, registry, state, data, preferredSelectedId) {
  var cards = validateCards(data.cards || [], registry, editorCollectionOptions(data));
  var selectedId = preferredSelectedId || state.selectedCardId;
  var displayCards;

  state.cards = cloneCards(cards);
  state.changeMetadataById = isObject(data.changeMetadataById) ? data.changeMetadataById : {};
  state.gitStatus = isObject(data.gitStatus) ? data.gitStatus : null;
  state.evidenceStatus = isObject(data.evidenceStatus) ? data.evidenceStatus : null;
  state.eligibilityById = isObject(data.eligibilityById) ? data.eligibilityById : {};
  state.changedFiles = Array.isArray(data.changedFiles) ? data.changedFiles : [];
  state.collectionVersion = typeof data.collectionVersion === "string" ? data.collectionVersion : null;
  state.diagnostics = Array.isArray(data.diagnostics) ? data.diagnostics : [];
  state.writesAllowed = data.writesAllowed !== false;
  renderEvidenceStatusBanner(root, state);
  renderCollectionDiagnostics(root, state);
  renderGitStatusWarning(root, state);
  renderChangedFiles(root, state);
  displayCards = displayCardsForState(state, registry);
  state.navigationStack = state.navigationStack.filter(function (cardId) {
    return Boolean(cardById(state.cards, cardId));
  });

  if (!cardById(state.cards, selectedId)) {
    selectedId = null;
  }
  state.selectedCardId = selectedId;
  state.currentPage = selectedId ? pageForCard(displayCards, selectedId, state.cardsPerPage) : 1;

  renderCardTablePage(root, registry, state, state.errorMessage);
  updateSelectedCardView(root, registry, state);
  updateEditAction(root, state);
  if (state.activeModal === "edit-card") {
    updateEditModal(root, registry, state);
  }
}

function setDeleteError(root, message) {
  var errorElement = root.querySelector("[data-delete-error]");

  if (errorElement) {
    errorElement.textContent = message || "";
    errorElement.hidden = !message;
  }
}

function updateEditAction(root, state) {
  var editButton = root.querySelector("[data-edit-card-trigger]");
  var deleteButton = root.querySelector("[data-delete-card]");
  var selectedCard = cardById(state.cards, state.selectedCardId);
  var newButton = root.querySelector('[data-open-modal="new-card"]');

  if (!selectedCard) {
    state.selectedCardId = null;
  }

  if (editButton) {
    editButton.disabled = !selectedCard || state.writesAllowed === false;
  }
  if (deleteButton) {
    deleteButton.disabled = !selectedCard || state.writesAllowed === false;
  }
  if (newButton) {
    newButton.disabled = state.writesAllowed === false;
  }
}

function updateEditModal(root, registry, state) {
  var body = root.querySelector("[data-edit-card-body]");
  var selectedCard = cardById(state.cards, state.selectedCardId);

  if (!body) {
    return;
  }

  if (!selectedCard) {
    body.innerHTML = '<div class="cm-modal-placeholder">Select a card to edit it.</div>';
    return;
  }

  body.innerHTML = renderEditorForm(registry, selectedCard, "edit", state);
  bindCardPickers(root, registry, state);
  bindSearchDropdowns(body);
  bindProductionBasisControls(body);
}

function updateBackAction(root, state) {
  var backButton = root.querySelector("[data-card-back]");
  var depth = root.querySelector("[data-nav-depth-slot]");

  if (backButton) {
    backButton.disabled = !state.navigationStack.length;
  }
  if (depth) {
    depth.innerHTML = renderNavigationDepth(state.navigationStack.length);
    depth.hidden = !state.navigationStack.length;
  }
}

function bindCardDetailNavigation(root, registry, state) {
  root.querySelectorAll("[data-card-link]").forEach(function (link) {
    link.addEventListener("click", function (event) {
      if (plainLeftClick(event)) {
        event.preventDefault();
        navigateToCard(root, registry, state, link.getAttribute("data-card-link"), { pushHistory: true });
      }
    });
  });
}

function updateSelectedCardView(root, registry, state) {
  var selectedCard = cardById(state.cards, state.selectedCardId);
  var view = root.querySelector("[data-card-view]");

  if (!selectedCard) {
    state.selectedCardId = null;
  }

  if (view) {
    view.innerHTML = renderSelectedCard(selectedCard, registry, state);
    bindCardDetailNavigation(root, registry, state);
  }

  updateBackAction(root, state);
}

function setSelectedCard(root, registry, state, selectedId) {
  var selectedCard = cardById(state.cards, selectedId);
  var rows = root.querySelectorAll(".cm-card-row");

  state.selectedCardId = selectedCard ? selectedId : null;
  updateCardUrl(state.selectedCardId);

  rows.forEach(function (row) {
    row.classList.toggle("is-selected", row.getAttribute("data-card-id") === selectedId);
  });

  updateSelectedCardView(root, registry, state);
  updateEditAction(root, state);
}

function selectCard(root, registry, state, selectedId) {
  var displayCards;

  if (!cardById(state.cards, selectedId)) {
    return false;
  }

  displayCards = displayCardsForState(state, registry);
  state.selectedCardId = selectedId;
  updateCardUrl(state.selectedCardId);
  state.currentPage = pageForCard(displayCards, selectedId, state.cardsPerPage);
  renderCardTablePage(root, registry, state, state.errorMessage);
  updateSelectedCardView(root, registry, state);
  updateEditAction(root, state);
  return true;
}

function navigateToCard(root, registry, state, selectedId, options) {
  var settings = options || {};

  if (!cardById(state.cards, selectedId)) {
    return false;
  }

  if (settings.clearHistory) {
    state.navigationStack = [];
  }
  if (settings.pushHistory && state.selectedCardId && state.selectedCardId !== selectedId) {
    state.navigationStack.push(state.selectedCardId);
  }

  return selectCard(root, registry, state, selectedId);
}

function navigateBack(root, registry, state) {
  var previousId;

  while (state.navigationStack.length) {
    previousId = state.navigationStack.pop();
    if (cardById(state.cards, previousId)) {
      selectCard(root, registry, state, previousId);
      return true;
    }
  }

  updateBackAction(root, state);
  return false;
}

function resetNewCardForm(root, registry, state) {
  var body = root.querySelector('[data-modal="new-card"] .cm-modal-body');

  if (body) {
    body.innerHTML = renderEditorForm(registry, null, "new", state);
    bindCardPickers(root, registry, state);
    bindSearchDropdowns(body);
    bindProductionBasisControls(body);
  }
}

function commitEditorForm(root, registry, state, mode) {
  var form = root.querySelector('[data-card-form="' + mode + '"]');
  var selectedCard = cardById(state.cards, state.selectedCardId);
  var commitButton = root.querySelector('[data-modal-commit="' + mode + '"]');
  var card;
  var existingCard;
  var validationCards;
  var method;
  var path;

  if (!state.isInternal || !form || state.isSaving || state.writesAllowed === false) {
    return false;
  }

  setFormError(form, "");
  try {
    card = collectEditorCard(form, registry, mode === "edit" ? selectedCard : null);
    validationCards = state.cards.filter(function (existing) {
      return existing.id !== card.id;
    }).concat([card]);
    validateCardForEditor(card, registry, validationCards, { checkReciprocal: false });

    if (mode === "new") {
      existingCard = cardById(state.cards, card.id);
      if (existingCard) {
        throw new Error("A card with ID '" + card.id + "' already exists.");
      }
    } else {
      if (!selectedCard) {
        throw new Error("Select a card before saving changes.");
      }
      if (card.id !== selectedCard.id) {
        throw new Error("Card ID cannot be changed while editing.");
      }
    }

  } catch (error) {
    setFormError(form, error.message);
    return false;
  }

  method = mode === "new" ? "POST" : "PUT";
  path = mode === "new" ? "/api/cifa-capabilities/cards" : "/api/cifa-capabilities/cards/" + encodeURIComponent(card.id);
  state.isSaving = true;
  if (commitButton) {
    commitButton.disabled = true;
    commitButton.textContent = mode === "new" ? "Creating..." : "Saving...";
  }

  sendEditorRequest(method, path, { card: card, collectionVersion: state.collectionVersion }).then(function (data) {
    applyEditorState(root, registry, state, data, card.id);
    setActiveModal(root, registry, state, null);
    if (mode === "new") {
      resetNewCardForm(root, registry, state);
    }
  }).catch(function (error) {
    setFormError(form, error.message || "Could not save card.");
  }).then(function () {
    state.isSaving = false;
    if (commitButton) {
      commitButton.disabled = false;
      commitButton.textContent = mode === "new" ? "Create" : "Save";
    }
  });

  return true;
}

function deleteSelectedCard(root, registry, state) {
  var selectedCard = cardById(state.cards, state.selectedCardId);
  var deleteButton = root.querySelector("[data-confirm-delete-card]");

  if (!state.isInternal || !selectedCard || state.isDeleting || state.writesAllowed === false) {
    return;
  }

  state.isDeleting = true;
  setDeleteError(root, "");
  if (deleteButton) {
    deleteButton.disabled = true;
    deleteButton.textContent = "Deleting...";
  }

  sendEditorRequest("DELETE", "/api/cifa-capabilities/cards/" + encodeURIComponent(selectedCard.id), { collectionVersion: state.collectionVersion }).then(function (data) {
    applyEditorState(root, registry, state, data, null);
    setActiveModal(root, registry, state, null);
  }).catch(function (error) {
    setDeleteError(root, error.message || "Could not delete card.");
  }).then(function () {
    state.isDeleting = false;
    if (deleteButton) {
      deleteButton.disabled = false;
      deleteButton.textContent = "Delete";
    }
  });
}

function setActiveModal(root, registry, state, modalName) {
  var frame = root.querySelector(".cm-app-frame");
  var backdrop = root.querySelector("[data-modal-backdrop]");
  var modalOpen = Boolean(modalName);

  state.activeModal = modalName;

  if (frame) {
    frame.classList.toggle("is-modal-open", modalOpen);
  }
  if (backdrop) {
    backdrop.setAttribute("aria-hidden", modalOpen ? "false" : "true");
  }

  root.querySelectorAll("[data-modal]").forEach(function (modal) {
    var isOpen = modal.getAttribute("data-modal") === modalName;
    modal.classList.toggle("is-open", isOpen);
    modal.setAttribute("aria-hidden", isOpen ? "false" : "true");
  });

  root.querySelectorAll("[data-open-modal]").forEach(function (button) {
    var isExpanded = button.getAttribute("data-open-modal") === modalName;
    button.setAttribute("aria-expanded", isExpanded ? "true" : "false");
  });

  if (modalName === "edit-card") {
    updateEditModal(root, registry, state);
  }
  if (modalName === "new-card") {
    resetNewCardForm(root, registry, state);
  }
  if (modalName === "delete-card") {
    setDeleteError(root, "");
  }
}

function syncSelectedRows(root, state) {
  root.querySelectorAll(".cm-card-row").forEach(function (row) {
    row.classList.toggle("is-selected", row.getAttribute("data-card-id") === state.selectedCardId);
  });
}

function bindCardRows(root, registry, state) {
  root.querySelectorAll(".cm-card-row").forEach(function (row) {
    row.addEventListener("click", function (event) {
      var cardId = row.getAttribute("data-card-id");

      if (plainLeftClick(event)) {
        event.preventDefault();
        state.navigationStack = [];
        setSelectedCard(root, registry, state, cardId);
      }
    });

    if (state.isInternal) {
      row.addEventListener("dblclick", function (event) {
        var cardId = row.getAttribute("data-card-id");

        event.preventDefault();
        state.navigationStack = [];
        setSelectedCard(root, registry, state, cardId);
        if (cardById(state.cards, cardId)) {
          setActiveModal(root, registry, state, "edit-card");
        }
      });
    }
  });
}

function bindPaginationControls(root, registry, state, errorMessage) {
  root.querySelectorAll("[data-card-page]").forEach(function (button) {
    button.addEventListener("click", function () {
      if (!button.disabled) {
        state.currentPage = parseInt(button.getAttribute("data-card-page"), 10);
        renderCardTablePage(root, registry, state, errorMessage);
      }
    });
  });
}

function updateCardsPerPageControls(root, state) {
  root.querySelectorAll("[data-cards-per-page-display]").forEach(function (input) {
    input.value = state.cardsPerPage;
  });
  root.querySelectorAll("[data-cards-per-page-step]").forEach(function (button) {
    var step = parseInt(button.getAttribute("data-cards-per-page-step"), 10);
    var nextValue = state.cardsPerPage + step;

    button.disabled = nextValue < MIN_CARDS_PER_PAGE;
  });
}

function setCardsPerPage(root, registry, state, cardsPerPage) {
  var displayCards;

  state.cardsPerPage = normalizeCardsPerPage(cardsPerPage);
  saveStoredCardsPerPagePreference(state.cardsPerPage);
  updateCardsPerPageControls(root, state);

  if (cardById(state.cards, state.selectedCardId)) {
    displayCards = displayCardsForState(state, registry);
    state.currentPage = pageForCard(displayCards, state.selectedCardId, state.cardsPerPage);
  }

  renderCardTablePage(root, registry, state, state.errorMessage);
}

function bindCardsPerPageControl(root, registry, state) {
  root.querySelectorAll("[data-cards-per-page-step]").forEach(function (button) {
    button.addEventListener("click", function () {
      var step = parseInt(button.getAttribute("data-cards-per-page-step"), 10);

      if (!button.disabled) {
        setCardsPerPage(root, registry, state, state.cardsPerPage + step);
      }
    });
  });
  updateCardsPerPageControls(root, state);
}

function updateOrderByChangesFirstControls(root, state) {
  root.querySelectorAll("[data-order-by-changes-first]").forEach(function (input) {
    input.checked = Boolean(state.orderByChangesFirst);
  });
}

function setOrderByChangesFirst(root, registry, state, enabled) {
  var displayCards;

  state.orderByChangesFirst = Boolean(enabled);
  saveStoredOrderByChangesFirstPreference(state.orderByChangesFirst);
  updateOrderByChangesFirstControls(root, state);

  if (cardById(state.cards, state.selectedCardId)) {
    displayCards = displayCardsForState(state, registry);
    state.currentPage = pageForCard(displayCards, state.selectedCardId, state.cardsPerPage);
  } else {
    state.currentPage = 1;
  }

  renderCardTablePage(root, registry, state, state.errorMessage);
}

function bindOrderByChangesFirstControl(root, registry, state) {
  root.querySelectorAll("[data-order-by-changes-first]").forEach(function (input) {
    input.addEventListener("change", function () {
      setOrderByChangesFirst(root, registry, state, input.checked);
    });
  });
  updateOrderByChangesFirstControls(root, state);
}

function updateCardSearchControls(root, state) {
  var filterMode = Boolean(state.filterPanelActive && hasActiveCardFilters(state));

  root.querySelectorAll("[data-card-search-input]").forEach(function (input) {
    var displayValue = filterMode ? "" : state.searchQuery;
    var placeholder = filterMode ? "Clear Filters" : "Search";

    if (input.value !== displayValue) {
      input.value = displayValue;
    }
    input.disabled = filterMode;
    input.setAttribute("placeholder", placeholder);
    input.setAttribute("aria-label", filterMode ? "Clear filters with the adjacent button" : "Search capability titles");
  });
  root.querySelectorAll("[data-filter-capability]").forEach(function (input) {
    if (input.value !== state.searchQuery) {
      input.value = state.searchQuery;
    }
  });
  root.querySelectorAll("[data-card-search-clear]").forEach(function (button) {
    button.hidden = !filterMode && !state.searchQuery;
    button.setAttribute("aria-label", filterMode ? "Clear filters" : "Clear search");
  });
  updateFilterControlsFromState(root, state);
}

function updateFilterControlsFromState(root, state) {
  FILTER_FIELD_KEYS.forEach(function (key) {
    var valueInput = root.querySelector('[name="filter__' + key + '"]');
    var dropdown = valueInput && valueInput.closest("[data-search-dropdown]");
    var textInput = dropdown && dropdown.querySelector("[data-search-dropdown-input]");
    var value = state.filters && state.filters[key] ? state.filters[key] : "";

    if (valueInput && valueInput.value !== value) {
      valueInput.value = value;
    }
    if (textInput && textInput.value !== value) {
      textInput.value = value;
    }
    if (dropdown) {
      updateSearchDropdownClear(dropdown);
    }
  });
}

function setCardSearchQuery(root, registry, state, query) {
  var displayCards;

  state.filters = createEmptyFilters();
  state.filterPanelActive = false;
  state.searchQuery = String(query || "");
  updateCardSearchControls(root, state);

  displayCards = displayCardsForState(state, registry);
  if (cardById(displayCards, state.selectedCardId)) {
    state.currentPage = pageForCard(displayCards, state.selectedCardId, state.cardsPerPage);
  } else {
    state.currentPage = 1;
  }

  renderCardTablePage(root, registry, state, state.errorMessage);
}

function clearAllFilters(root, registry, state) {
  var displayCards;

  state.searchQuery = "";
  state.filters = createEmptyFilters();
  state.filterPanelActive = false;
  updateCardSearchControls(root, state);

  displayCards = displayCardsForState(state, registry);
  if (cardById(displayCards, state.selectedCardId)) {
    state.currentPage = pageForCard(displayCards, state.selectedCardId, state.cardsPerPage);
  } else {
    state.currentPage = 1;
  }

  renderCardTablePage(root, registry, state, state.errorMessage);
}

function bindCardSearchControl(root, registry, state) {
  root.querySelectorAll("[data-card-search-input]").forEach(function (input) {
    input.addEventListener("input", function () {
      setCardSearchQuery(root, registry, state, input.value);
    });
  });
  root.querySelectorAll("[data-card-search-clear]").forEach(function (button) {
    button.addEventListener("click", function () {
      var input = root.querySelector("[data-card-search-input]");

      clearAllFilters(root, registry, state);
      if (input) {
        input.focus();
      }
    });
  });
  updateCardSearchControls(root, state);
}

function filterControlValue(root, key) {
  var control = root.querySelector('[name="filter__' + key + '"]');

  return control ? control.value.trim() : "";
}

function applyFilterControls(root, registry, state) {
  var capabilityControl = root.querySelector("[data-filter-capability]");
  var displayCards;

  state.searchQuery = capabilityControl ? capabilityControl.value.trim() : "";
  FILTER_FIELD_KEYS.forEach(function (key) {
    state.filters[key] = filterControlValue(root, key);
  });
  state.filterPanelActive = hasActiveCardFilters(state);
  updateCardSearchControls(root, state);

  displayCards = displayCardsForState(state, registry);
  if (cardById(displayCards, state.selectedCardId)) {
    state.currentPage = pageForCard(displayCards, state.selectedCardId, state.cardsPerPage);
  } else {
    state.currentPage = 1;
  }

  renderCardTablePage(root, registry, state, state.errorMessage);
  setActiveModal(root, registry, state, null);
}

function bindFilterControls(root, registry, state) {
  root.querySelectorAll("[data-apply-filters]").forEach(function (button) {
    button.addEventListener("click", function () {
      applyFilterControls(root, registry, state);
    });
  });
}

function renderCapabilityFilter(state) {
  return [
    '<label class="cm-field cm-filter-field">',
    '  <span>Capability</span>',
    '  <input type="search" placeholder="Search" value="' + escapeHtml(state.searchQuery || "") + '" data-filter-capability autocomplete="off">',
    '</label>'
  ].join("");
}

function renderEnumFilter(registry, key, state) {
  var field = registry.fieldByKey[key];

  if (!field || !Array.isArray(field.options)) {
    return "";
  }

  return [
    '<label class="cm-field cm-filter-field">',
    '  <span>' + escapeHtml(field.label) + '</span>',
       renderSearchDropdown("filter__" + key, field.options, state.filters[key] || "", "All"),
    '</label>'
  ].join("");
}

function renderFilterControls(registry, state) {
  return [
    renderCapabilityFilter(state),
    renderEnumFilter(registry, "tool", state),
    renderEnumFilter(registry, "status", state),
    renderEnumFilter(registry, "minAbaqusVersion", state),
    renderEnumFilter(registry, "minCoreformIgaVersion", state)
  ].filter(Boolean).join("");
}

function renderOrderByChangesFirstFilter(isInternal, state) {
  if (!isInternal) {
    return "";
  }

  return [
    '    <div class="cm-filter-bottom-options">',
    '      <label class="cm-filter-checkbox">',
    '        <input type="checkbox" data-order-by-changes-first' + (state.orderByChangesFirst ? " checked" : "") + '>',
    '        <span>Order by changes first</span>',
    '      </label>',
    '    </div>'
  ].join("");
}

function renderCardSearchControl(state) {
  var filterMode = Boolean(state.filterPanelActive && hasActiveCardFilters(state));
  var query = filterMode ? "" : state.searchQuery || "";
  var placeholder = filterMode ? "Clear Filters" : "Search";

  return [
    '<label class="cm-card-search">',
    '  <span class="cm-sr-only">Search capability titles</span>',
    '  <input type="search" placeholder="' + escapeHtml(placeholder) + '" value="' + escapeHtml(query) + '" data-card-search-input autocomplete="off" aria-label="' + (filterMode ? "Clear filters with the adjacent button" : "Search capability titles") + '"' + (filterMode ? " disabled" : "") + '>',
    '  <button class="cm-search-clear cm-card-search-clear" type="button" data-card-search-clear aria-label="' + (filterMode ? "Clear filters" : "Clear search") + '"' + (filterMode || query ? "" : " hidden") + '>×</button>',
    '</label>'
  ].join("");
}

function renderCardTablePage(root, registry, state, errorMessage) {
  var list = root.querySelector("[data-card-list]");
  var footer = root.querySelector("[data-card-pagination]");
  var displayCards = displayCardsForState(state, registry);
  var pagination = paginateCards(displayCards, state.currentPage, state.cardsPerPage);

  state.currentPage = pagination.currentPage;

  if (list) {
    list.innerHTML = renderCardRows(pagination.visibleCards, registry, errorMessage, state.changeMetadataById, hasActiveCardFilters(state));
  }
  if (footer) {
    footer.innerHTML = renderPagination(pagination, errorMessage);
  }

  bindCardRows(root, registry, state);
  bindPaginationControls(root, registry, state, errorMessage);
  syncSelectedRows(root, state);
  updateEditAction(root, state);
}

function render(root, data, errorMessage) {
  var mode = root.getAttribute("data-mode") || "public";
  var isInternal = mode === "internal";
  var appTitle = root.getAttribute("data-title") || "Coreform IGA for Abaqus Capabilities";
  var cards = data && Array.isArray(data.cards) ? data.cards : [];
  var registry = data && data.registry ? data.registry : createEmptyRegistry();
  var displayCards;
  var initialCardId;
  var state = {
    activeModal: null,
    cards: cloneCards(cards),
    cardsPerPage: loadStoredCardsPerPagePreference(),
    changeMetadataById: data && isObject(data.changeMetadataById) ? data.changeMetadataById : {},
    currentPage: 1,
    errorMessage: errorMessage,
    evidenceStatus: data && isObject(data.evidenceStatus) ? data.evidenceStatus : null,
    eligibilityById: data && isObject(data.eligibilityById) ? data.eligibilityById : {},
    changedFiles: data && Array.isArray(data.changedFiles) ? data.changedFiles : [],
    collectionVersion: data && typeof data.collectionVersion === "string" ? data.collectionVersion : null,
    gitStatus: data && isObject(data.gitStatus) ? data.gitStatus : null,
    diagnostics: data && Array.isArray(data.diagnostics) ? data.diagnostics : [],
    writesAllowed: !data || data.writesAllowed !== false,
    filterPanelActive: false,
    filters: createEmptyFilters(),
    isInternal: isInternal,
    isDeleting: false,
    isSaving: false,
    navigationStack: [],
    orderByChangesFirst: isInternal ? loadStoredOrderByChangesFirstPreference() : false,
    searchQuery: "",
    selectedCardId: null
  };

  displayCards = displayCardsForState(state, registry);
  initialCardId = cardIdFromLocation();
  state.selectedCardId = cardById(state.cards, initialCardId) ? initialCardId : null;
  if (state.selectedCardId) {
    state.currentPage = pageForCard(displayCards, state.selectedCardId, state.cardsPerPage);
  }

  root.innerHTML = [
    '<div class="cm-app-frame">',
    '  <div class="cm-page-backdrop" data-modal-backdrop aria-hidden="true"></div>',
    '  <section class="cm-modal-panel cm-filter-panel" data-modal="filter" aria-hidden="true">',
    '    <div class="cm-filter-header">',
    '      <h2>Filters</h2>',
    '      <p>Filter the Capability List by capability text and registry field values.</p>',
    '    </div>',
    '    <div class="cm-filter-grid">',
           renderFilterControls(registry, state),
    '    </div>',
         renderOrderByChangesFirstFilter(isInternal, state),
    '    <div class="cm-filter-footer">',
    '      <label class="cm-cards-per-page-control">',
    '        <span>Cards per page</span>',
    '        <span class="cm-cards-per-page-stepper">',
    '          <input type="text" value="' + state.cardsPerPage + '" data-cards-per-page-display disabled>',
    '          <button class="cm-cards-per-page-button" type="button" data-cards-per-page-step="-' + CARDS_PER_PAGE_STEP + '" aria-label="Show fewer cards per page">-</button>',
    '          <button class="cm-cards-per-page-button" type="button" data-cards-per-page-step="' + CARDS_PER_PAGE_STEP + '" aria-label="Show more cards per page">+</button>',
    '        </span>',
    '      </label>',
    '      <div class="cm-filter-actions">',
    '        <button class="cm-button cm-button-primary" type="button" data-apply-filters>Apply</button>',
    '        <button class="cm-button cm-button-secondary" type="button" data-modal-close>Close</button>',
    '      </div>',
    '    </div>',
    '  </section>',
       renderHelpModal(isInternal),
       renderInternalModals(isInternal, registry),
    '  <section class="cm-titlebar" aria-label="Page title">',
    '    <a class="cm-brand-link" href="' + COREFORM_HOME_URL + '" aria-label="Coreform home">',
    '      <img class="cm-brand-logo" src="' + COREFORM_LOGO_URL + '" width="333" height="96" alt="Coreform">',
    '    </a>',
    '    <div class="cm-page-title-lockup">',
    '      <img class="cm-title-logo" src="' + TITLE_LOGO_URL + '" width="67" height="67" alt="">',
    '      <h1>' + escapeHtml(appTitle) + '</h1>',
    '    </div>',
    '    <div class="cm-titlebar-actions">',
           renderHelpAction(),
    '    </div>',
    '  </section>',
    '  <section class="cm-evidence-status" data-evidence-status hidden aria-live="polite"></section>',
    '  <section class="cm-evidence-status cm-evidence-status--invalid" data-collection-diagnostics hidden aria-live="polite"></section>',
    '  <section class="cm-evidence-status" data-changed-files hidden aria-live="polite"></section>',
    '  <section class="cm-evidence-status cm-evidence-status--invalid" data-git-status-warning hidden aria-live="polite"></section>',
    '  <main class="cm-workspace">',
    '    <section class="cm-panel cm-card-table-panel" aria-label="Capability List">',
    '      <div class="cm-panel-header">',
    '        <h1>Capability List</h1>',
    '        <div class="cm-list-header-actions">',
               renderCardSearchControl(state),
    '          <button class="cm-filter-toggle" type="button" data-open-modal="filter" data-filter-toggle aria-expanded="false">Filters</button>',
    '        </div>',
    '      </div>',
           renderCardTableHeader(registry),
    '      <div class="cm-card-list" data-card-list></div>',
    '      <footer class="cm-card-pagination" data-card-pagination aria-live="polite"></footer>',
    '    </section>',
    '    <section class="cm-panel cm-card-view-panel" aria-label="Capability View">',
    '      <div class="cm-panel-header">',
    '        <h1>Capability View</h1>',
    '        <div class="cm-card-view-actions">',
               renderNewCardAction(isInternal),
               renderEditCardAction(isInternal),
    '          <div class="cm-card-view-navigation">',
    '            <span data-nav-depth-slot></span>',
    '            <button class="cm-button cm-button-secondary cm-card-back-button" type="button" data-card-back disabled>Back</button>',
    '          </div>',
    '        </div>',
    '      </div>',
    '      <div class="cm-card-view" data-card-view>',
             renderSelectedCard(cardById(state.cards, state.selectedCardId), registry, state),
    '      </div>',
    '    </section>',
    '  </main>',
    '</div>'
  ].join("");

  renderEvidenceStatusBanner(root, state);
  renderCollectionDiagnostics(root, state);
  renderGitStatusWarning(root, state);
  bindCardPickers(root, registry, state);
  bindSearchDropdowns(root);
  bindProductionBasisControls(root);
  bindGlobalDismissHandlers();
  bindCardsPerPageControl(root, registry, state);
  bindCardSearchControl(root, registry, state);
  bindFilterControls(root, registry, state);
  bindOrderByChangesFirstControl(root, registry, state);

  root.querySelectorAll("[data-open-modal]").forEach(function (button) {
    button.addEventListener("click", function () {
      if (!button.disabled) {
        setActiveModal(root, registry, state, button.getAttribute("data-open-modal"));
      }
    });
  });

  root.querySelectorAll("[data-modal-close]").forEach(function (button) {
    button.addEventListener("click", function () {
      setActiveModal(root, registry, state, null);
    });
  });

  root.querySelectorAll("[data-modal-commit]").forEach(function (button) {
    button.addEventListener("click", function () {
      commitEditorForm(root, registry, state, button.getAttribute("data-modal-commit"));
    });
  });

  root.querySelector("[data-modal-backdrop]").addEventListener("click", function () {
    setActiveModal(root, registry, state, null);
  });

  root.querySelectorAll("[data-card-back]").forEach(function (button) {
    button.addEventListener("click", function () {
      if (!button.disabled) {
        navigateBack(root, registry, state);
      }
    });
  });

  root.querySelectorAll("[data-delete-card]").forEach(function (button) {
    button.addEventListener("click", function () {
      if (!button.disabled) {
        setActiveModal(root, registry, state, "delete-card");
      }
    });
  });

  root.querySelectorAll("[data-confirm-delete-card]").forEach(function (button) {
    button.addEventListener("click", function () {
      deleteSelectedCard(root, registry, state);
    });
  });

  renderCardTablePage(root, registry, state, errorMessage);
  updateSelectedCardView(root, registry, state);
}

function init() {
  var root = document.getElementById("cifa-capabilities-app");
  var mode = root ? root.getAttribute("data-mode") || "public" : "public";

  if (!root || root.getAttribute("data-cifa-capabilities-initialized") === "true") {
    return false;
  }

  root.setAttribute("data-cifa-capabilities-initialized", "true");
  if (root) {
    loadData(mode).then(function (data) {
      render(root, data, null);
    }).catch(function (error) {
      render(root, { cards: [], registry: createEmptyRegistry(), changeMetadataById: {} }, error.message);
    });
  }

  return true;
}

function scheduleInit() {
  var attempts = 0;
  var observer = null;
  var timer;

  if (init()) {
    return;
  }

  function stopWatching() {
    if (timer) {
      window.clearInterval(timer);
    }
    if (observer) {
      observer.disconnect();
    }
  }

  if (window.MutationObserver && document.documentElement) {
    observer = new MutationObserver(function () {
      if (init()) {
        stopWatching();
      }
    });
    observer.observe(document.documentElement, {
      childList: true,
      subtree: true
    });
  }

  timer = window.setInterval(function () {
    attempts += 1;
    if (init() || attempts > 4800) {
      stopWatching();
    }
  }, 25);

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", function () {
      if (init()) {
        stopWatching();
      }
    });
  }
}

scheduleInit();
