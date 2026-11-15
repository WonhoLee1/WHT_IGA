// Form collection, searchable dropdowns, and relationship picker behavior.

function setFormError(form, message) {
  var errorElement = form && form.querySelector("[data-form-error]");

  if (!errorElement) {
    if (message) {
      window.alert(message);
    }
    return;
  }

  errorElement.textContent = message || "";
  errorElement.hidden = !message;
}

function selectedCardPickerIds(picker) {
  return selectedCardPickerEntries(picker).map(function (entry) {
    return entry.id;
  });
}

function selectedCardPickerIdsInForm(picker) {
  var form = picker.closest("form");
  var selected = Object.create(null);

  if (!form) {
    selectedCardPickerIds(picker).forEach(function (cardId) {
      selected[cardId] = true;
    });
    return selected;
  }

  form.querySelectorAll("[data-card-picker-field]").forEach(function (formPicker) {
    selectedCardPickerIds(formPicker).forEach(function (cardId) {
      selected[cardId] = true;
    });
  });

  return selected;
}

function selectedCardPickerEntries(picker) {
  return Array.prototype.map.call(picker.querySelectorAll("[data-card-picker-selected]"), function (row) {
    return {
      id: row.getAttribute("data-card-picker-selected"),
      supported: row.getAttribute("data-card-picker-supported") !== "false"
    };
  });
}

function cardPickerCurrentId(picker) {
  var form = picker.closest("form");
  var idInput = form && form.elements.id;

  return idInput ? idInput.value.trim() : "";
}

function renderCardPickerMenu(picker, registry, state) {
  var input = picker.querySelector("[data-card-picker-input]");
  var menu = picker.querySelector("[data-card-picker-menu]");
  var query = input ? input.value.trim().toLowerCase() : "";
  var selected = Object.create(null);
  var currentId = cardPickerCurrentId(picker);
  var matches;

  if (!menu) {
    return;
  }
  updateCardPickerSearchClear(picker);

  selected = selectedCardPickerIdsInForm(picker);

  matches = state.cards.filter(function (card) {
    var haystack = (card.id + " " + titleForCard(card, registry)).toLowerCase();

    return card.id !== currentId && !selected[card.id] && (!query || haystack.indexOf(query) !== -1);
  });

  if (!matches.length) {
    menu.innerHTML = '<div class="cm-card-picker-empty cm-search-dropdown-empty">No matching cards.</div>';
    menu.hidden = false;
    return;
  }

  menu.innerHTML = matches.map(function (card) {
    return [
      '<button class="cm-card-picker-option cm-search-dropdown-option" type="button" data-card-picker-add="' + escapeHtml(card.id) + '">',
      '  <span class="cm-card-picker-option-id">' + escapeHtml(card.id) + '</span>',
      '  <span class="cm-card-picker-option-title">' + escapeHtml(titleForCard(card, registry)) + '</span>',
      '</button>'
    ].join("");
  }).join("");
  menu.hidden = false;
}

function renderSelectedCardPickerRows(picker, registry, state) {
  var list = picker.querySelector("[data-card-picker-selected-list]");
  var entries;

  if (!list) {
    return;
  }

  entries = selectedCardPickerEntries(picker);
  list.innerHTML = renderCardPickerEditorRows(entries, registry, state);
}

function addCardPickerId(picker, cardId, registry, state) {
  var list = picker.querySelector("[data-card-picker-selected-list]");
  var input = picker.querySelector("[data-card-picker-input]");
  var entries = selectedCardPickerEntries(picker);

  if (!list || selectedCardPickerIds(picker).indexOf(cardId) !== -1) {
    return;
  }

  entries.push({
    id: cardId,
    supported: true
  });
  list.innerHTML = renderCardPickerEditorRows(entries, registry, state);
  if (input) {
    input.value = "";
  }
  updateCardPickerSearchClear(picker);
}

function updateCardPickerSearchClear(picker) {
  var input = picker && picker.querySelector("[data-card-picker-input]");
  var clearButton = picker && picker.querySelector("[data-card-picker-clear]");

  if (clearButton) {
    clearButton.hidden = !input || !input.value.trim();
  }
}

function closeCardPickerMenu(picker) {
  var menu = picker && picker.querySelector("[data-card-picker-menu]");

  if (menu) {
    menu.hidden = true;
  }
}

function closeCardPickerMenus(container) {
  container.querySelectorAll("[data-card-picker-menu]").forEach(function (menu) {
    menu.hidden = true;
  });
}

function blurCardPickerInput(picker) {
  var input = picker && picker.querySelector("[data-card-picker-input]");

  if (input) {
    input.blur();
  }
  if (picker && document.activeElement && picker.contains(document.activeElement)) {
    document.activeElement.blur();
  }
}

function updateCurrentParentAction(picker, state) {
  var button = picker.querySelector("[data-card-picker-current-parent]");
  var selectedCard = state.selectedCardId ? cardById(state.cards, state.selectedCardId) : null;
  var alreadySelected = selectedCard && selectedCardPickerIdsInForm(picker)[selectedCard.id] === true;

  if (!button) {
    return;
  }

  button.disabled = !selectedCard || alreadySelected;
  button.textContent = "Make child of current card";
  if (!selectedCard) {
    button.title = "Select a card in Capability View first.";
  } else if (alreadySelected) {
    button.title = "The current card is already selected as a relationship.";
  } else {
    button.title = "";
  }
}

function updateCurrentParentActions(container, state) {
  container.querySelectorAll("[data-card-picker-field]").forEach(function (picker) {
    updateCurrentParentAction(picker, state);
  });
}

function bindCardPickers(root, registry, state) {
  root.querySelectorAll("[data-card-picker-field]").forEach(function (picker) {
    var input = picker.querySelector("[data-card-picker-input]");
    var menu = picker.querySelector("[data-card-picker-menu]");
    var form = picker.closest("form");
    var idInput = form && form.elements.id;

    if (picker.getAttribute("data-card-picker-bound") === "true") {
      return;
    }
    picker.setAttribute("data-card-picker-bound", "true");

    if (input) {
      input.addEventListener("focus", function () {
        renderCardPickerMenu(picker, registry, state);
      });
      input.addEventListener("input", function () {
        renderCardPickerMenu(picker, registry, state);
      });
    }

    if (idInput) {
      idInput.addEventListener("input", function () {
        if (menu && !menu.hidden) {
          renderCardPickerMenu(picker, registry, state);
        }
      });
    }

    if (form && form.getAttribute("data-card-picker-dismiss-bound") !== "true") {
      form.setAttribute("data-card-picker-dismiss-bound", "true");
      form.addEventListener("click", function (event) {
        if (!event.target.closest("[data-card-picker-input]")) {
          closeCardPickerMenus(form);
        }
      });
    }

    picker.addEventListener("click", function (event) {
      var addButton = event.target.closest("[data-card-picker-add]");
      var clearButton = event.target.closest("[data-card-picker-clear]");
      var removeButton = event.target.closest("[data-card-picker-remove]");
      var supportButton = event.target.closest("[data-card-picker-support-toggle]");
      var currentParentButton = event.target.closest("[data-card-picker-current-parent]");
      var supportRow;

      if (clearButton && input) {
        input.value = "";
        input.focus();
        renderCardPickerMenu(picker, registry, state);
        return;
      }
      if (addButton) {
        addCardPickerId(picker, addButton.getAttribute("data-card-picker-add"), registry, state);
        renderCardPickerMenu(picker, registry, state);
        closeCardPickerMenu(picker);
        blurCardPickerInput(picker);
      }
      if (removeButton) {
        removeButton.closest("[data-card-picker-selected]").remove();
        renderSelectedCardPickerRows(picker, registry, state);
        closeCardPickerMenu(picker);
      }
      if (supportButton) {
        supportRow = supportButton.closest("[data-card-picker-selected]");
        supportRow.setAttribute("data-card-picker-supported", supportRow.getAttribute("data-card-picker-supported") === "false" ? "true" : "false");
        renderSelectedCardPickerRows(picker, registry, state);
        closeCardPickerMenu(picker);
      }
      if (currentParentButton && state.selectedCardId) {
        addCardPickerId(picker, state.selectedCardId, registry, state);
        renderCardPickerMenu(picker, registry, state);
        closeCardPickerMenu(picker);
        blurCardPickerInput(picker);
      }
      updateCurrentParentActions(form || picker, state);
    });

    updateCurrentParentActions(form || picker, state);
    updateCardPickerSearchClear(picker);
  });
}

function searchDropdownMatches(option, query) {
  return !query || option.toLowerCase().indexOf(query) !== -1;
}

function closeSearchDropdown(dropdown) {
  var menu = dropdown && dropdown.querySelector("[data-search-dropdown-menu]");

  if (menu) {
    menu.hidden = true;
  }
}

function updateSearchDropdownValue(dropdown) {
  var input = dropdown.querySelector("[data-search-dropdown-input]");
  var valueInput = dropdown.querySelector("[data-search-dropdown-value]");
  var query = input ? input.value.trim() : "";
  var exactValue = "";

  dropdown.querySelectorAll("[data-search-dropdown-option]").forEach(function (optionButton) {
    var option = optionButton.getAttribute("data-search-dropdown-option") || "";

    if (option.toLowerCase() === query.toLowerCase()) {
      exactValue = option;
    }
  });

  if (valueInput) {
    valueInput.value = exactValue;
  }
  updateSearchDropdownClear(dropdown);
}

function normalizeSearchDropdownText(dropdown) {
  var input = dropdown.querySelector("[data-search-dropdown-input]");
  var valueInput = dropdown.querySelector("[data-search-dropdown-value]");

  if (!input || !valueInput) {
    return;
  }

  if (valueInput.value) {
    input.value = valueInput.value;
  } else {
    input.value = "";
  }
  updateSearchDropdownClear(dropdown);
}

function updateSearchDropdownClear(dropdown) {
  var input = dropdown.querySelector("[data-search-dropdown-input]");
  var valueInput = dropdown.querySelector("[data-search-dropdown-value]");
  var clearButton = dropdown.querySelector("[data-search-dropdown-clear]");
  var hasValue = Boolean((valueInput && valueInput.value) || (input && input.value.trim()));

  if (clearButton) {
    clearButton.hidden = !hasValue;
  }
}

function renderSearchDropdownMenu(dropdown) {
  var input = dropdown.querySelector("[data-search-dropdown-input]");
  var menu = dropdown.querySelector("[data-search-dropdown-menu]");
  var empty = dropdown.querySelector("[data-search-dropdown-empty]");
  var query = input ? input.value.trim().toLowerCase() : "";
  var visibleCount = 0;

  if (!menu) {
    return;
  }

  dropdown.querySelectorAll("[data-search-dropdown-option]").forEach(function (optionButton) {
    var option = optionButton.getAttribute("data-search-dropdown-option") || "";
    var visible = searchDropdownMatches(option, query);

    optionButton.hidden = !visible;
    if (visible) {
      visibleCount += 1;
    }
  });

  if (empty) {
    empty.hidden = visibleCount !== 0;
  }
  updateSearchDropdownClear(dropdown);
  menu.hidden = false;
}

function setSearchDropdownSelection(dropdown, value) {
  var input = dropdown.querySelector("[data-search-dropdown-input]");
  var valueInput = dropdown.querySelector("[data-search-dropdown-value]");

  if (valueInput) {
    valueInput.value = value;
  }
  if (input) {
    input.value = value;
    input.blur();
  }
  updateSearchDropdownClear(dropdown);
  closeSearchDropdown(dropdown);
}

function clearSearchDropdownSelection(dropdown) {
  var input = dropdown.querySelector("[data-search-dropdown-input]");
  var valueInput = dropdown.querySelector("[data-search-dropdown-value]");

  if (valueInput) {
    valueInput.value = "";
  }
  if (input) {
    input.value = "";
    input.focus();
  }
  updateSearchDropdownClear(dropdown);
  renderSearchDropdownMenu(dropdown);
}

function selectSearchDropdownText(input) {
  window.setTimeout(function () {
    if (document.activeElement === input) {
      input.select();
    }
  }, 0);
}

function bindSearchDropdowns(root) {
  root.querySelectorAll("[data-search-dropdown]").forEach(function (dropdown) {
    var input = dropdown.querySelector("[data-search-dropdown-input]");

    if (dropdown.getAttribute("data-search-dropdown-bound") === "true") {
      return;
    }
    dropdown.setAttribute("data-search-dropdown-bound", "true");

    if (input) {
      input.addEventListener("focus", function () {
        renderSearchDropdownMenu(dropdown);
        selectSearchDropdownText(input);
      });
      input.addEventListener("input", function () {
        updateSearchDropdownValue(dropdown);
        renderSearchDropdownMenu(dropdown);
      });
      input.addEventListener("blur", function () {
        window.setTimeout(function () {
          if (!dropdown.contains(document.activeElement)) {
            normalizeSearchDropdownText(dropdown);
          }
        }, 0);
      });
    }

    dropdown.addEventListener("click", function (event) {
      var optionButton = event.target.closest("[data-search-dropdown-option]");
      var clearButton = event.target.closest("[data-search-dropdown-clear]");

      if (clearButton) {
        clearSearchDropdownSelection(dropdown);
        return;
      }
      if (optionButton) {
        setSearchDropdownSelection(dropdown, optionButton.getAttribute("data-search-dropdown-option") || "");
      }
    });

    updateSearchDropdownClear(dropdown);
  });
}

function bindGlobalDismissHandlers() {
  if (globalDismissHandlersBound) {
    return;
  }
  globalDismissHandlersBound = true;

  document.addEventListener("click", function (event) {
    document.querySelectorAll("[data-card-picker-field]").forEach(function (picker) {
      if (!picker.contains(event.target)) {
        closeCardPickerMenu(picker);
      }
    });
    document.querySelectorAll("[data-search-dropdown]").forEach(function (dropdown) {
      if (!dropdown.contains(event.target)) {
        closeSearchDropdown(dropdown);
      }
    });
  });
}

function nonEmptyLines(value) {
  return String(value || "").split(/\r?\n/).map(function (line) {
    return line.trim();
  }).filter(Boolean);
}

function bindProductionBasisControls(root) {
  root.querySelectorAll("[data-production-basis]").forEach(function (fieldset) {
    var select = fieldset.querySelector('select[name="productionBasis__kind"]');
    if (select.dataset.basisBound) {
      return;
    }
    select.dataset.basisBound = "true";
    select.addEventListener("change", function () {
      fieldset.querySelectorAll("[data-basis-panel]").forEach(function (panel) {
        panel.hidden = panel.getAttribute("data-basis-panel") !== select.value;
      });
    });
  });
}

function collectEditorCard(form, registry, seedCard) {
  var card = seedCard ? cloneCard(seedCard) : {};

  registry.fields.forEach(function (field) {
    if (field.generated === true) {
      delete card[field.key];
    }
  });

  registry.editorFields.forEach(function (field) {
    var control = form.elements[field.key];
    var picker;
    var selectedEntries;

    if (field.editor.widget === "production-basis") {
      var kind = form.elements.productionBasis__kind.value;
      delete card.productionBasis;
      delete card.evidenceContract;
      if (kind === "verification") {
        card.productionBasis = { kind: kind };
        try {
          card.evidenceContract = JSON.parse(form.elements.productionBasis__contract.value);
        } catch (error) {
          throw new Error("Verification contract must be valid JSON: " + error.message);
        }
      } else if (kind === "override") {
        card.productionBasis = { kind: kind, justification: form.elements.productionBasis__justification.value.trim() };
        var references = nonEmptyLines(form.elements.productionBasis__references.value);
        if (references.length) {
          card.productionBasis.references = references;
        }
      }
    } else if (field.type === "card-list") {
      picker = form.querySelector('[data-card-picker-field="' + field.key + '"]');
      selectedEntries = picker ? selectedCardPickerEntries(picker) : [];
      if (selectedEntries.length) {
        card[field.key] = selectedEntries;
      } else {
        delete card[field.key];
      }
    } else if (field.type === "string-list" || field.type === "url-list") {
      if (control && nonEmptyLines(control.value).length) {
        card[field.key] = nonEmptyLines(control.value);
      } else {
        delete card[field.key];
      }
    } else if (field.type === "plan") {
      var stageControl = form.elements[field.key + "__stage"];
      var monthControl = form.elements[field.key + "__month"];
      var yearControl = form.elements[field.key + "__year"];
      var stage = stageControl ? stageControl.value.trim() : "";
      var month = monthControl ? monthControl.value.trim() : "";
      var year = yearControl ? yearControl.value.trim() : "";
      var date = month && year ? month + " " + year : "";

      if (stage || month || year) {
        card[field.key] = { stage: stage, date: date };
      } else {
        delete card[field.key];
      }
    } else if (field.type === "boolean-map") {
      var anyBlocked = false;
      var blockedValues = {};

      field.items.forEach(function (item) {
        var itemControl = form.elements[field.key + "__" + item.key];
        var checked = Boolean(itemControl && itemControl.checked);

        blockedValues[item.key] = checked;
        anyBlocked = anyBlocked || checked;
      });
      if (anyBlocked) {
        card[field.key] = blockedValues;
      } else {
        delete card[field.key];
      }
    } else if (field.type === "grouped-string-lists") {
      var groupedValues = {};
      var hasGroupedValues = false;

      field.groups.forEach(function (group) {
        var groupControl = form.elements[field.key + "__" + group.key];
        var lines = groupControl ? nonEmptyLines(groupControl.value) : [];

        if (lines.length) {
          groupedValues[group.key] = lines;
          hasGroupedValues = true;
        }
      });
      if (hasGroupedValues) {
        card[field.key] = groupedValues;
      } else {
        delete card[field.key];
      }
    } else if (control) {
      if (!field.required && control.value.trim() === "") {
        delete card[field.key];
      } else {
        card[field.key] = control.value.trim();
      }
    }
  });

  return card;
}
