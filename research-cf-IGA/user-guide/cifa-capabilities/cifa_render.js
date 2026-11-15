// Markup renderers for the capability list, capability view, help, and editor modals.

function evidenceStatusMessage(status) {
  var producing = String(status.producingRevision || "").slice(0, 12);
  var current = String(status.currentRevision || "").slice(0, 12);
  var required = String(status.requiredRevision || status.currentRevision || "").slice(0, 12);
  var override = required !== current ? "; required source revision " + required + " (CIFA_SOURCE_COMMIT); checkout HEAD " + current : "";
  var generated = status.generatedAt ? ", generated " + status.generatedAt : "";

  if (status.state === "current-exact") {
    return "Current evidence — index from revision " + producing + generated + override + ".";
  }
  if (status.state === "stale-unverified") {
    return "Stale (unverified) evidence — index from revision " + producing + generated +
      "; this checkout is at " + current + override +
      ". " + (status.detail || "Evidence-backed eligibility is unverified until the index is regenerated.");
  }
  if (status.state === "missing") {
    return "Evidence pending — no capability evidence index was found. Structural editing is unaffected.";
  }
  var message = "Evidence index invalid — " + (status.detail || "the index was rejected.");
  if (status.issues && status.issues.length) {
    message += " First issue: " + status.issues[0];
  }
  return message;
}

function renderEvidenceStatusBanner(root, state) {
  var element = root.querySelector("[data-evidence-status]");

  if (!element) {
    return;
  }
  if (!state.evidenceStatus || !state.evidenceStatus.state) {
    element.hidden = true;
    return;
  }
  element.className = "cm-evidence-status cm-evidence-status--" + state.evidenceStatus.state;
  element.textContent = evidenceStatusMessage(state.evidenceStatus);
  element.hidden = false;
}

function renderGitStatusWarning(root, state) {
  var element = root.querySelector("[data-git-status-warning]");
  if (!element) {
    return;
  }
  var unavailable = state.gitStatus && state.gitStatus.state === "unavailable";
  element.hidden = !unavailable;
  element.textContent = unavailable ? state.gitStatus.detail + " Structural editing remains available; missing change badges do not mean a clean worktree." : "";
}

function statusPillClass(value) {
  switch (String(value || "").toLowerCase()) {
    case "production":
      return "cm-status-pill-production";
    case "beta":
      return "cm-status-pill-beta";
    case "alpha":
      return "cm-status-pill-alpha";
    case "unsupported":
      return "cm-status-pill-unsupported";
    default:
      return "";
  }
}

function renderCollectionDiagnostics(root, state) {
  var element = root.querySelector("[data-collection-diagnostics]");
  if (!element) {
    return;
  }
  var diagnostics = state.diagnostics || [];
  element.hidden = !diagnostics.length;
  if (!diagnostics.length) {
    element.textContent = "";
    return;
  }
  element.innerHTML = '<strong>Source repair needed — Create, Save, and Delete are paused.</strong>' +
    '<p>Valid cards remain available for inspection. Fix these files under data/cards/ in your text editor, then reload. No source was changed automatically.</p>' +
    '<ul>' + diagnostics.map(function (item) {
      return '<li><code>' + escapeHtml(item.file) + '</code>: ' + escapeHtml(item.message) + '</li>';
    }).join("") + '</ul><button type="button" class="cm-button" data-reload-sources>Reload sources</button>';
  element.querySelector("[data-reload-sources]").addEventListener("click", function () { window.location.reload(); });
}

function renderStatusPill(value) {
  var classes = "cm-status-pill";
  var statusClass = statusPillClass(value);

  if (statusClass) {
    classes += " " + statusClass;
  }

  return '<span class="' + classes + '">' + escapeHtml(value) + '</span>';
}

function renderTableCell(card, field, changeMetadataById) {
  var value = fieldValue(card, field);
  var classes = "cm-card-cell";
  var style = field.table && field.table.style;
  var isTitle = field.detail && field.detail.role === "title";
  var badge = isTitle ? changeMetadataBadge(changeMetadataForCard(changeMetadataById, card.id)) : "";

  if (isTitle) {
    classes += " cm-card-cell-main";
  }

  if (style === "badge") {
    return '<span class="' + classes + '">' + renderStatusPill(value) + '</span>';
  }

  if (badge) {
    return [
      '<span class="' + classes + ' cm-card-cell-with-badge">',
           badge,
      '  <span class="cm-card-cell-text">' + escapeHtml(value) + '</span>',
      '</span>'
    ].join("");
  }

  return '<span class="' + classes + '">' + escapeHtml(value) + '</span>';
}

function renderCardTableHeader(registry) {
  if (!registry.tableFields.length) {
    return "";
  }

  return [
    '<div class="cm-card-table-header" aria-hidden="true"' + tableGridStyle(registry) + '>',
         registry.tableFields.map(function (field) {
           return '<span>' + escapeHtml(field.label) + '</span>';
         }).join(""),
    '</div>'
  ].join("");
}

function renderCardRows(cards, registry, errorMessage, changeMetadataById, hasActiveFilters) {
  if (errorMessage) {
    return [
      '<div class="cm-card-data-message cm-card-data-error">',
      '  <strong>Card data could not load.</strong>',
      '  <span>' + escapeHtml(errorMessage) + '</span>',
      '</div>'
    ].join("");
  }

  if (!cards.length) {
    return '<div class="cm-card-data-message">' + (hasActiveFilters ? "No matching cards." : "No cards available.") + '</div>';
  }

  return cards.map(function (card) {
    return [
      '<a class="cm-card-row" href="' + escapeHtml(cardHref(card.id)) + '" data-card-id="' + escapeHtml(card.id) + '"' + tableGridStyle(registry) + '>',
           registry.tableFields.map(function (field) {
             return renderTableCell(card, field, changeMetadataById);
           }).join(""),
      '</a>'
    ].join("");
  }).join("");
}

function totalCardPages(cards, cardsPerPage) {
  return Math.max(1, Math.ceil(cards.length / normalizeCardsPerPage(cardsPerPage)));
}

function clampPage(page, totalPages) {
  var numericPage = parseInt(page, 10);

  if (!numericPage || numericPage < 1) {
    return 1;
  }
  if (numericPage > totalPages) {
    return totalPages;
  }

  return numericPage;
}

function paginateCards(cards, page, cardsPerPage) {
  var totalCards = cards.length;
  var pageSize = normalizeCardsPerPage(cardsPerPage);
  var totalPages = totalCardPages(cards, pageSize);
  var currentPage = clampPage(page, totalPages);
  var startIndex = totalCards ? (currentPage - 1) * pageSize : 0;
  var endIndex = Math.min(startIndex + pageSize, totalCards);

  return {
    currentPage: currentPage,
    endIndex: endIndex,
    totalCards: totalCards,
    totalPages: totalPages,
    startIndex: startIndex,
    visibleCards: cards.slice(startIndex, endIndex)
  };
}

function renderPaginationButton(label, page, disabled, current, ariaLabel) {
  var classes = "cm-pagination-button" + (current ? " is-current" : "");

  return [
    '<button class="' + classes + '" type="button" data-card-page="' + escapeHtml(page) + '"',
         disabled ? " disabled" : "",
         current ? ' aria-current="page"' : "",
         ariaLabel ? ' aria-label="' + escapeHtml(ariaLabel) + '"' : "",
    '>',
         escapeHtml(label),
    '</button>'
  ].join("");
}

function renderPagination(pagination, errorMessage) {
  var buttons = [];
  var page;

  if (errorMessage || !pagination.totalCards) {
    return "";
  }

  buttons.push(renderPaginationButton("First", 1, pagination.currentPage === 1, false));
  buttons.push(renderPaginationButton("←", pagination.currentPage - 1, pagination.currentPage === 1, false, "Previous page"));

  for (page = 1; page <= pagination.totalPages; page += 1) {
    buttons.push(renderPaginationButton(page, page, false, page === pagination.currentPage));
  }

  buttons.push(renderPaginationButton("→", pagination.currentPage + 1, pagination.currentPage === pagination.totalPages, false, "Next page"));
  buttons.push(renderPaginationButton("Last", pagination.totalPages, pagination.currentPage === pagination.totalPages, false));

  return [
    '<div class="cm-pagination-info">',
    '  Showing ' + (pagination.startIndex + 1) + '-' + pagination.endIndex + ' of ' + pagination.totalCards,
    '</div>',
    '<div class="cm-pagination-buttons" aria-label="Card table pages">',
         buttons.join(""),
    '</div>'
  ].join("");
}

function titleForCard(card, registry) {
  var titleField = fieldByDetailRole(registry, "title") || registry.fieldByKey.id;

  return displayFieldValue(card, titleField);
}

function renderTextSegment(value) {
  return escapeHtml(value).replace(/\n/g, "<br>");
}

function parsedTextLink(source, state, registry) {
  var parsed = parsedInlineLinkTarget(source);
  var targetCard;

  if (!parsed) {
    return null;
  }

  if (parsed.type === "external") {
    return {
      type: "external",
      url: parsed.url,
      label: parsed.label,
      found: true
    };
  }

  targetCard = state ? cardById(state.cards, parsed.id) : null;
  return {
    type: "card",
    id: parsed.id,
    label: parsed.label || (targetCard ? titleForCard(targetCard, registry) : parsed.id),
    found: Boolean(targetCard)
  };
}

function renderTextWithCardLinks(value, state, registry) {
  var text = String(value || "");
  var output = [];
  var lastIndex = 0;

  text.replace(/\[\[([^\]]+)\]\]/g, function (match, source, offset) {
    var link;

    output.push(renderTextSegment(text.slice(lastIndex, offset)));
    link = parsedTextLink(source, state, registry);
    if (link && link.found) {
      if (link.type === "external") {
        output.push(
          '<a class="cm-card-link cm-url-link" href="' + escapeHtml(link.url) + '" target="_blank" rel="noopener noreferrer">' +
            escapeHtml(link.label) +
          '</a>'
        );
      } else {
        output.push(
          '<a class="cm-card-link" href="' + escapeHtml(cardHref(link.id)) + '" data-card-link="' + escapeHtml(link.id) + '">' +
            escapeHtml(link.label) +
          '</a>'
        );
      }
    } else {
      output.push('<span class="cm-card-link-missing">' + escapeHtml(match) + '</span>');
    }
    lastIndex = offset + match.length;
    return match;
  });

  output.push(renderTextSegment(text.slice(lastIndex)));
  return output.join("");
}

function renderSummarySection(card, field, state, registry) {
  if (!field) {
    return "";
  }

  return [
    '<section class="cm-detail-summary-section">',
    '  <h3>' + escapeHtml(field.label) + '</h3>',
    '  <p class="cm-detail-summary">' + renderTextWithCardLinks(fieldValue(card, field), state, registry) + '</p>',
    '</section>'
  ].join("");
}

function appendUniqueRelationship(entries, seen, entry, currentId) {
  if (entry.id && entry.id !== currentId && !seen[entry.id]) {
    seen[entry.id] = true;
    entries.push({
      id: entry.id,
      supported: entry.supported === true
    });
  }
}

function sortRelationshipEntries(entries, state, registry) {
  return entries.slice().sort(function (left, right) {
    var leftCard = state ? cardById(state.cards, left.id) : null;
    var rightCard = state ? cardById(state.cards, right.id) : null;
    var leftTitle = leftCard ? titleForCard(leftCard, registry) : left.id;
    var rightTitle = rightCard ? titleForCard(rightCard, registry) : right.id;

    if (left.supported !== right.supported) {
      return left.supported ? -1 : 1;
    }
    return leftTitle.localeCompare(rightTitle) || left.id.localeCompare(right.id);
  });
}

function relatedRelationships(card, cards, relationKey, inverseKey) {
  var entries = [];
  var seen = Object.create(null);

  cardListValue(card, relationKey).forEach(function (entry) {
    appendUniqueRelationship(entries, seen, entry, card.id);
  });

  cards.forEach(function (otherCard) {
    var inverseEntry = relationshipEntryFor(cardListValue(otherCard, inverseKey), card.id);

    if (inverseEntry) {
      appendUniqueRelationship(entries, seen, { id: otherCard.id, supported: inverseEntry.supported }, card.id);
    }
  });

  return entries;
}

function renderSupportedMark(supported) {
  return [
    '<span class="cm-supported-mark ' + (supported ? "cm-supported-mark-yes" : "cm-supported-mark-no") + '"',
    ' aria-label="' + (supported ? "Supported" : "Not supported") + '">',
         supported ? "✓" : "×",
    '</span>'
  ].join("");
}

function renderRelatedRow(entry, state, registry) {
  var card = cardById(state.cards, entry.id);
  var toolField = registry.fieldByKey.tool;
  var statusField = registry.fieldByKey.status;

  if (!card) {
    return "";
  }
  return [
    '<tr>',
    '  <td>' + renderSupportedMark(entry.supported) + '</td>',
    '  <td><a class="cm-card-link" href="' + escapeHtml(cardHref(card.id)) + '" data-card-link="' + escapeHtml(card.id) + '">' + escapeHtml(titleForCard(card, registry)) + '</a></td>',
    '  <td>' + (toolField ? escapeHtml(displayFieldValue(card, toolField)) : "") + '</td>',
    '  <td>' + (statusField ? renderStatusPill(displayFieldValue(card, statusField)) : "") + '</td>',
    '</tr>'
  ].join("");
}

function renderRelatedGroup(label, entries, state, registry) {
  var sortedEntries = sortRelationshipEntries(entries, state, registry).filter(function (entry) {
    return Boolean(cardById(state.cards, entry.id));
  });

  if (!sortedEntries.length) {
    return "";
  }

  return [
    '<section class="cm-related-group">',
    '  <h4>' + escapeHtml(label) + '</h4>',
    '  <div class="cm-relationship-table-wrap">',
    '    <table class="cm-relationship-table cm-related-table">',
    '      <thead><tr><th>Supported</th><th>Capability</th><th>Tool</th><th>Status</th></tr></thead>',
    '      <tbody>',
             sortedEntries.map(function (entry) {
               return renderRelatedRow(entry, state, registry);
             }).join(""),
    '      </tbody>',
    '    </table>',
    '  </div>',
    '</section>'
  ].join("");
}

function renderRelatedCards(card, state, registry) {
  var children = relatedRelationships(card, state.cards, "children", "parents");
  var parents = relatedRelationships(card, state.cards, "parents", "children");

  if (!children.length && !parents.length) {
    return "";
  }

  return [
    '<section class="cm-related-section">',
    '  <h3>Related Cards</h3>',
       renderRelatedGroup("Children", children, state, registry),
       renderRelatedGroup("Parents", parents, state, registry),
    '</section>'
  ].join("");
}

function renderSimpleList(items, className, renderer) {
  var values = Array.isArray(items) ? items : [];

  if (!values.length) {
    return "";
  }

  return [
    '<ul class="' + className + '">',
         values.map(function (item) {
           return '<li>' + renderer(item) + '</li>';
         }).join(""),
    '</ul>'
  ].join("");
}

function renderBooleanMap(value, field) {
  var map = isObject(value) ? value : {};

  return [
    '<div class="cm-boolean-map">',
         field.items.map(function (item) {
           var checked = map[item.key] === true;

           return [
             '<span class="cm-boolean-map-item">',
             '  <input type="checkbox" aria-disabled="true" tabindex="-1"' + (checked ? " checked" : "") + '>',
             '  <span>' + escapeHtml(item.label) + '</span>',
             '</span>'
           ].join("");
         }).join(""),
    '</div>'
  ].join("");
}

function renderGroupedStringLists(value, field) {
  var groups = isObject(value) ? value : {};
  var output = field.groups.map(function (group) {
    var items = Array.isArray(groups[group.key]) ? groups[group.key] : [];

    if (!items.length) {
      return "";
    }

    return [
      '<div class="cm-grouped-list">',
      '  <strong>' + escapeHtml(group.label) + '</strong>',
         renderSimpleList(items, "cm-string-list", escapeHtml),
      '</div>'
    ].join("");
  }).filter(Boolean);

  return output.length ? '<div class="cm-grouped-lists">' + output.join("") + '</div>' : "";
}

function renderDetailFieldValue(card, field) {
  var value = fieldValue(card, field);

  if (value === undefined || value === "" || value === null) {
    return "";
  }
  if (field.type === "url-list") {
    return renderSimpleList(value, "cm-url-list", function (url) {
      return '<a href="' + escapeHtml(url) + '" target="_blank" rel="noopener noreferrer">' + escapeHtml(url) + '</a>';
    });
  }
  if (field.type === "string-list") {
    return renderSimpleList(value, "cm-string-list", escapeHtml);
  }
  if (field.type === "plan") {
    return isObject(value) ? escapeHtml([value.stage, value.date].filter(Boolean).join(" - ")) : "";
  }
  if (field.type === "boolean-map") {
    return renderBooleanMap(value, field);
  }
  if (field.type === "grouped-string-lists") {
    return renderGroupedStringLists(value, field);
  }

  return escapeHtml(displayFieldValue(card, field));
}

function renderDetailList(card, fields) {
  var rows = fields.map(function (field) {
    var renderedValue = renderDetailFieldValue(card, field);

    return renderedValue ? '<div><dt>' + escapeHtml(field.label) + '</dt><dd>' + renderedValue + '</dd></div>' : "";
  }).join("");

  return rows ? '<dl class="cm-detail-list">' + rows + '</dl>' : "";
}

function renderInternalOnlyDetails(card, fields) {
  var detailList = renderDetailList(card, fields);

  if (!detailList) {
    return "";
  }

  return [
    '<section class="cm-internal-detail-section">',
    '  <h3>Internal Only</h3>',
       detailList,
    '</section>'
  ].join("");
}

function renderInternalCardId(card, registry, state) {
  var idField = registry.fieldByKey.id;
  var cardId = idField ? displayFieldValue(card, idField) : "";

  if (!state || !state.isInternal || !cardId) {
    return "";
  }

  return [
    '<div class="cm-detail-id">',
    '  <span>ID</span>',
    '  <code>' + escapeHtml(cardId) + '</code>',
    '</div>'
  ].join("");
}

function renderCardEligibility(card, state) {
  var eligibility = state && state.isInternal && state.eligibilityById && state.eligibilityById[card.id];
  var labels = { "schema-valid": "Schema-valid", "evidence-pending": "Evidence pending", "publishable": "Publishable Production basis" };
  if (!eligibility) {
    return "";
  }
  return '<section class="cm-card-eligibility" data-card-eligibility="' + escapeHtml(eligibility.state) + '">' +
    '<h3>Saved card: ' + escapeHtml(labels[eligibility.state] || eligibility.state) + '</h3>' +
    '<p>' + escapeHtml(eligibility.detail) + '</p><p class="cm-form-hint">' + escapeHtml(eligibility.scope) + '</p></section>';
}

function renderChangedFiles(root, state) {
  var element = root.querySelector("[data-changed-files]");
  if (!element) {
    return;
  }
  var files = state.isInternal && state.changedFiles || [];
  element.hidden = !files.length;
  element.innerHTML = files.length ? '<h2>Saved locally — changed files</h2><ul>' + files.map(function (file) {
    return '<li><code>' + escapeHtml(file) + '</code></li>';
  }).join("") + '</ul><p>Review the Git diff and commit these changes. Saving does not run tests or publish.</p>' : "";
}

function renderSelectedCard(card, registry, state) {
  var titleField;
  var badgeField;
  var summaryField;
  var metadataFields;
  var publicMetadataFields;
  var internalMetadataFields;

  if (!card) {
    return '<div class="cm-empty-view">Click on a Card to View</div>';
  }

  titleField = fieldByDetailRole(registry, "title") || registry.fieldByKey.id;
  badgeField = fieldByDetailRole(registry, "badge");
  summaryField = fieldByDetailRole(registry, "summary");
  metadataFields = registry.detailFields.filter(function (field) {
    return (field.detail.role || "metadata") === "metadata";
  });
  publicMetadataFields = metadataFields.filter(function (field) {
    return field.internalOnly !== true;
  });
  internalMetadataFields = metadataFields.filter(function (field) {
    return field.internalOnly === true && field.key !== "id";
  });

  return [
    '<article class="cm-card-detail">',
    '  <h2 class="cm-detail-title">' + escapeHtml(displayFieldValue(card, titleField)) + '</h2>',
         renderInternalCardId(card, registry, state),
         badgeField ? '  ' + renderStatusPill(displayFieldValue(card, badgeField)) : "",
         renderSummarySection(card, summaryField, state, registry),
         renderDetailList(card, publicMetadataFields),
         state ? renderRelatedCards(card, state, registry) : "",
         state && state.isInternal ? renderInternalOnlyDetails(card, internalMetadataFields) : "",
         renderCardEligibility(card, state),
    '</article>'
  ].join("");
}

function renderNewCardAction(isInternal) {
  if (!isInternal) {
    return "";
  }

  return '<button class="cm-button cm-button-secondary cm-action-button" type="button" data-open-modal="new-card" aria-expanded="false">New Card</button>';
}

function renderEditCardAction(isInternal) {
  if (!isInternal) {
    return "";
  }

  return '<button class="cm-button cm-button-secondary cm-action-button" type="button" data-open-modal="edit-card" data-edit-card-trigger aria-expanded="false" disabled>Edit Card</button>';
}

function renderHelpAction() {
  return '<button class="cm-help-button" type="button" data-open-modal="help" aria-label="Help" aria-expanded="false">?</button>';
}

function renderPublicHelpContent() {
  return [
    '<div class="cm-help-body">',
    '  <section>',
    '    <h3>Cards</h3>',
    '    <p>Each row in the Capability List represents one capability card. Select a row to view the card details in the Capability View panel.</p>',
    '  </section>',
    '  <section>',
    '    <h3>Relationships</h3>',
    '    <p>Cards can point to related cards through the Children and Parents sections. Select a related card to navigate to it.</p>',
    '  </section>',
    '  <section>',
    '    <h3>Navigation</h3>',
    '    <p>Use the page buttons at the bottom of the Capability List to move through cards. Card links and related cards jump directly to the target card, external links open supporting websites, and the Back button returns to the previous card in that path.</p>',
    '  </section>',
    '  <section>',
    '    <h3>Filters</h3>',
    '    <p>The Filters panel narrows the Capability List by title text, tool, status, and minimum required versions. It also controls how many cards appear per page.</p>',
    '  </section>',
    '</div>'
  ].join("");
}

function renderInternalHelpContent() {
  return [
    '<div class="cm-help-body">',
    '  <section>',
    '    <h3>Top-right actions</h3>',
    '    <ul>',
    '      <li><strong>New Card</strong> writes a new card JSON file into the local repo served by the editor.</li>',
    '      <li><strong>Edit Card</strong> updates the selected local card file. Double-clicking a table row opens the same editor.</li>',
    '      <li><strong>Delete</strong> is available inside Edit Card and removes the selected card after confirmation.</li>',
    '    </ul>',
    '  </section>',
    '  <section>',
    '    <h3>Local changes</h3>',
    '    <p>Create, Save, and Delete write directly to the checked-out card JSON files on disk. Use Git in your terminal for commits, resets, branches, pushes, and merge requests.</p>',
    '  </section>',
    '  <section>',
    '    <h3>Links and relationships</h3>',
    '    <p>Use <code>[[card-id]]</code> to link to a card by title, or <code>[[label &lt;card-id&gt;]]</code> to choose the linked text.</p>',
    '    <p>Use <code>[[label &lt;https://example.com&gt;]]</code> to add an external summary link, or <code>[[https://example.com]]</code> to use the URL as the linked text.</p>',
    '    <p>Children and Parents relationships are mirrored automatically. Adding one side updates the other card so both JSON files stay in sync.</p>',
    '  </section>',
    '  <section>',
    '    <h3>Badges</h3>',
    '    <p><strong>New</strong> marks untracked or newly added card files. <strong>Edited</strong> marks modified tracked card files.</p>',
    '  </section>',
    '</div>'
  ].join("");
}

function renderHelpModal(isInternal) {
  return [
    '  <section class="cm-modal-panel cm-action-modal cm-help-modal" data-modal="help" aria-hidden="true">',
    '    <div class="cm-modal-header">',
    '      <h2>Help</h2>',
    '    </div>',
    '    <div class="cm-modal-body">',
           isInternal ? renderInternalHelpContent() : renderPublicHelpContent(),
    '    </div>',
    '    <div class="cm-modal-actions">',
    '      <button class="cm-button cm-button-secondary" type="button" data-modal-close>Close</button>',
    '    </div>',
    '  </section>'
  ].join("");
}

function renderNavigationDepth(stackDepth) {
  var visibleDepth = Math.min(stackDepth, NAV_DEPTH_DISPLAY_LIMIT);
  var overflow = stackDepth - visibleDepth;
  var chevrons = [];
  var index;

  if (!stackDepth) {
    return "";
  }

  for (index = 0; index < visibleDepth; index += 1) {
    chevrons.push(
      '<span class="cm-nav-chevron' + (index === 0 ? " cm-nav-chevron-active" : "") + '" aria-hidden="true"></span>'
    );
  }

  if (overflow > 0) {
    chevrons.push('<span class="cm-nav-depth-overflow" aria-hidden="true">+' + escapeHtml(overflow) + '</span>');
  }

  return [
    '<span class="cm-nav-depth" data-nav-depth aria-label="Navigation depth: ' + escapeHtml(stackDepth) + '">',
         chevrons.join(""),
    '</span>'
  ].join("");
}

function renderCardPickerEditorRow(entry, registry, state) {
  var card = state ? cardById(state.cards, entry.id) : null;
  var title = card ? titleForCard(card, registry) : entry.id;

  return [
    '<tr data-card-picker-selected="' + escapeHtml(entry.id) + '" data-card-picker-supported="' + escapeHtml(entry.supported ? "true" : "false") + '">',
    '  <td><button class="cm-supported-toggle ' + (entry.supported ? "is-supported" : "is-unsupported") + '" type="button" data-card-picker-support-toggle="' + escapeHtml(entry.id) + '" aria-label="' + (entry.supported ? "Mark not supported" : "Mark supported") + '">' + (entry.supported ? "✓" : "×") + '</button></td>',
    '  <td>' + escapeHtml(title) + '</td>',
    '  <td class="cm-relationship-id-cell">' + escapeHtml(entry.id) + '</td>',
    '  <td><button class="cm-card-picker-remove" type="button" data-card-picker-remove="' + escapeHtml(entry.id) + '" aria-label="Remove ' + escapeHtml(entry.id) + '">Remove</button></td>',
    '</tr>'
  ].join("");
}

function renderCardPickerEditorRows(entries, registry, state) {
  var sortedEntries = sortRelationshipEntries(entries, state, registry);

  if (!sortedEntries.length) {
    return '<tr class="cm-relationship-empty-row"><td colspan="4">No related cards selected.</td></tr>';
  }

  return sortedEntries.map(function (entry) {
    return renderCardPickerEditorRow(entry, registry, state);
  }).join("");
}

function renderCardPickerControl(field, card, mode, registry, state) {
  var values = cardListValue(card, field.key);
  var showCurrentParentAction = mode === "new" && field.key === "parents";

  return [
    '<div class="cm-form-field cm-card-picker" data-card-picker-field="' + escapeHtml(field.key) + '">',
    '  <span>' + escapeHtml(field.label) + '</span>',
    '  <div class="cm-card-picker-search cm-search-dropdown">',
    '    <input class="cm-search-dropdown-input" type="search" data-card-picker-input placeholder="Search card IDs" autocomplete="off">',
    '    <button class="cm-search-clear cm-card-picker-search-clear" type="button" data-card-picker-clear aria-label="Clear card search" hidden>×</button>',
    '    <div class="cm-card-picker-menu cm-search-dropdown-menu" data-card-picker-menu hidden></div>',
    '  </div>',
       showCurrentParentAction ? '<button class="cm-card-picker-current-parent" type="button" data-card-picker-current-parent>Make child of current card</button>' : "",
    '  <div class="cm-relationship-table-wrap cm-card-picker-table-wrap">',
    '    <table class="cm-relationship-table cm-card-picker-table">',
    '      <thead><tr><th>Supported</th><th>Capability</th><th>ID</th><th></th></tr></thead>',
    '      <tbody data-card-picker-selected-list>',
             renderCardPickerEditorRows(values, registry, state),
    '      </tbody>',
    '    </table>',
    '  </div>',
    '</div>'
  ].join("");
}

function listEditorValue(value) {
  return Array.isArray(value) ? value.join("\n") : "";
}

function renderSearchDropdown(name, options, value, placeholder) {
  var selectedValue = options.indexOf(value) === -1 ? "" : value;

  return [
    '<div class="cm-search-dropdown" data-search-dropdown>',
    '  <input type="hidden" name="' + escapeHtml(name) + '" value="' + escapeHtml(selectedValue) + '" data-search-dropdown-value>',
    '  <input class="cm-search-dropdown-input" type="search" value="' + escapeHtml(selectedValue) + '" placeholder="' + escapeHtml(placeholder || "Search options") + '" data-search-dropdown-input autocomplete="off" role="searchbox">',
    '  <button class="cm-search-clear cm-search-dropdown-clear" type="button" data-search-dropdown-clear aria-label="Clear selection"' + (selectedValue ? "" : " hidden") + '>×</button>',
    '  <div class="cm-search-dropdown-menu" data-search-dropdown-menu hidden role="listbox">',
         options.map(function (option) {
           return [
             '<button class="cm-search-dropdown-option" type="button" data-search-dropdown-option="' + escapeHtml(option) + '" role="option">',
             '  <span class="cm-search-dropdown-option-primary">' + escapeHtml(option) + '</span>',
             '</button>'
           ].join("");
         }).join(""),
    '    <div class="cm-search-dropdown-empty" data-search-dropdown-empty hidden>No matching options.</div>',
    '  </div>',
    '</div>'
  ].join("");
}

function renderPlanControl(field, card) {
  var value = isObject(fieldValue(card, field)) ? fieldValue(card, field) : {};
  var planDate = parsePlanDate(value.date);

  return [
    '<div class="cm-form-field cm-plan-field">',
    '  <span>' + escapeHtml(field.label) + '</span>',
    '  <div class="cm-plan-controls">',
         renderSearchDropdown(field.key + "__stage", field.options, value.stage || "", "Stage"),
         renderSearchDropdown(field.key + "__month", PLAN_MONTHS, planDate.month, "Month"),
         renderSearchDropdown(field.key + "__year", PLAN_YEARS, planDate.year, "Year"),
    '  </div>',
    '</div>'
  ].join("");
}

function renderBooleanMapControl(field, card) {
  var value = isObject(fieldValue(card, field)) ? fieldValue(card, field) : {};

  return [
    '<fieldset class="cm-form-field cm-boolean-map-editor">',
    '  <legend>' + escapeHtml(field.label) + '</legend>',
         field.items.map(function (item) {
           return [
             '<label>',
             '  <input type="checkbox" name="' + escapeHtml(field.key + "__" + item.key) + '"' + (value[item.key] === true ? " checked" : "") + '>',
             '  <span>' + escapeHtml(item.label) + '</span>',
             '</label>'
           ].join("");
         }).join(""),
    '</fieldset>'
  ].join("");
}

function renderGroupedStringListsControl(field, card) {
  var value = isObject(fieldValue(card, field)) ? fieldValue(card, field) : {};

  return [
    '<div class="cm-form-field cm-grouped-string-lists-editor">',
    '  <span>' + escapeHtml(field.label) + '</span>',
         field.groups.map(function (group) {
           return [
             '<label>',
             '  <span>' + escapeHtml(group.label) + '</span>',
             '  <textarea name="' + escapeHtml(field.key + "__" + group.key) + '" rows="4">' + escapeHtml(listEditorValue(value[group.key])) + '</textarea>',
             '</label>'
           ].join("");
         }).join(""),
    '</div>'
  ].join("");
}

function renderProductionBasisControl(card) {
  var basis = card.productionBasis || {};
  var kind = basis.kind || "";
  var example = { minimumLevel: "Exercised", anyOf: [{ name: "case-name", ctestPattern: "^Test_abaqus_verification_case_name$" }] };
  return [
    '<fieldset class="cm-production-basis" data-production-basis>',
    '<legend>Production basis</legend>',
    '<p class="cm-form-hint">Choose verification or an explicit override. A taxonomy parent with neither derives eligibility from all its supported children; every such child must be Production and eligible.</p>',
    '<label class="cm-form-field"><span>Basis</span><select name="productionBasis__kind">',
    [['', 'No independent basis / supported children'], ['verification', 'Verification Manual evidence'], ['override', 'Explicit Production override']].map(function (option) {
      return '<option value="' + option[0] + '"' + (kind === option[0] ? ' selected' : '') + '>' + option[1] + '</option>';
    }).join(""),
    '</select></label>',
    '<div data-basis-panel="verification"' + (kind === 'verification' ? '' : ' hidden') + '>',
    '<label class="cm-form-field"><span>Verification contract (JSON)</span>',
    '<textarea name="productionBasis__contract" rows="10" spellcheck="false" placeholder="' + escapeHtml(JSON.stringify(example, null, 2)) + '">' + escapeHtml(card.evidenceContract ? JSON.stringify(card.evidenceContract, null, 2) : '') + '</textarea></label>',
    '<p class="cm-form-hint">Use anyOf cases with a name and Python ctestPattern; optional contract-level minimumLevel is Exercised or Numerically validated. Existing nested matchers are preserved. Eligibility requires a passing official test represented in the Verification Manual. Save validates the contract but never runs it.</p>',
    '</div>',
    '<div data-basis-panel="override"' + (kind === 'override' ? '' : ' hidden') + '>',
    '<label class="cm-form-field"><span>Why this is Production</span><textarea name="productionBasis__justification" rows="4">' + escapeHtml(basis.justification || '') + '</textarea></label>',
    '<label class="cm-form-field"><span>References (optional, one HTTP(S) URL per line)</span><textarea name="productionBasis__references" rows="3">' + escapeHtml((basis.references || []).join('\n')) + '</textarea></label>',
    '<p class="cm-form-hint">A justification is sufficient: no test, reviewer, date, expiry, or attachment is required. This declaration is labeled as an override, never verification evidence.</p>',
    '</div>',
    '<p class="cm-form-hint">Switching basis replaces the previous declaration and contract on Save. Changing a basis does not change the status selector.</p>',
    '</fieldset>'
  ].join("");
}

function renderEditorControl(field, card, mode, registry, state) {
  var editor = field.editor || {};
  var value = fieldValue(card, field);
  var required = field.required ? " required" : "";
  var readOnly = mode === "edit" && editor.readOnlyOnEdit === true;
  var control;

  if (editor.widget === "production-basis") {
    return renderProductionBasisControl(card);
  }
  if (editor.widget === "card-picker") {
    return renderCardPickerControl(field, card, mode, registry, state);
  }
  if (editor.widget === "plan") {
    return renderPlanControl(field, card);
  }
  if (editor.widget === "boolean-map") {
    return renderBooleanMapControl(field, card);
  }
  if (editor.widget === "grouped-string-lists") {
    return renderGroupedStringListsControl(field, card);
  }

  if (editor.widget === "textarea") {
    control = [
      '<textarea name="' + escapeHtml(field.key) + '" rows="5"' + required + '>',
           escapeHtml(field.type === "string-list" || field.type === "url-list" ? listEditorValue(value) : value),
      '</textarea>'
    ].join("");
  } else if (editor.widget === "select") {
    control = renderSearchDropdown(field.key, field.options, value, field.label);
  } else {
    control = '<input type="text" name="' + escapeHtml(field.key) + '" value="' + escapeHtml(value) + '"' + required + (readOnly ? ' readonly aria-readonly="true"' : "") + '>';
  }

  return [
    '<label class="cm-form-field">',
    '  <span>' + escapeHtml(field.label) + '</span>',
       control,
    '</label>'
  ].join("");
}

function addEditorField(fields, seen, field) {
  if (field && field.editor && field.editor.visible !== false && !seen[field.key]) {
    seen[field.key] = true;
    fields.push(field);
  }
}

function orderedEditorFieldGroups(registry) {
  var seen = Object.create(null);
  var publicFields = [];
  var internalFields = [];
  var titleField = fieldByDetailRole(registry, "title");
  var badgeField = fieldByDetailRole(registry, "badge");
  var summaryField = fieldByDetailRole(registry, "summary");

  addEditorField(publicFields, seen, titleField);
  addEditorField(publicFields, seen, registry.fieldByKey.id);
  addEditorField(publicFields, seen, badgeField);
  addEditorField(publicFields, seen, summaryField);

  registry.detailFields.forEach(function (field) {
    var role = field.detail.role || "metadata";

    if (field.internalOnly !== true && role === "metadata") {
      addEditorField(publicFields, seen, field);
    }
  });
  registry.detailFields.forEach(function (field) {
    var role = field.detail.role || "metadata";

    if (field.internalOnly !== true && role === "relationship") {
      addEditorField(publicFields, seen, field);
    }
  });
  registry.editorFields.forEach(function (field) {
    if (field.internalOnly !== true) {
      addEditorField(publicFields, seen, field);
    }
  });
  registry.editorFields.forEach(function (field) {
    if (field.internalOnly === true && field.key !== "id") {
      addEditorField(internalFields, seen, field);
    }
  });

  return {
    publicFields: publicFields,
    internalFields: internalFields
  };
}

function renderEditorControls(fields, card, mode, registry, state) {
  return fields.map(function (field) {
    return renderEditorControl(field, card, mode, registry, state);
  }).join("");
}

function renderInternalEditorSection(fields, card, mode, registry, state) {
  var controls = renderEditorControls(fields, card, mode, registry, state);

  if (!controls) {
    return "";
  }

  return [
    '<section class="cm-form-section cm-form-section-internal">',
    '  <h3>Internal Only</h3>',
       controls,
    '</section>'
  ].join("");
}

function renderEditorForm(registry, card, mode, state) {
  var editorCard = card || {};
  var groups = orderedEditorFieldGroups(registry);

  if (!registry.editorFields.length) {
    return '<div class="cm-modal-placeholder">No editable fields are defined.</div>';
  }

  return [
    '<form class="cm-card-form" data-card-form="' + escapeHtml(mode) + '">',
    '  <div class="cm-form-error" data-form-error hidden></div>',
    '  <p class="cm-form-hint">Create/Save checks structure and writes source files locally. Eligibility shown outside this form describes the last saved card, not these unsaved edits.</p>',
       renderEditorControls(groups.publicFields, editorCard, mode, registry, state),
       renderInternalEditorSection(groups.internalFields, editorCard, mode, registry, state),
    '</form>'
  ].join("");
}

function renderInternalModals(isInternal, registry) {
  if (!isInternal) {
    return "";
  }

  return [
    '  <section class="cm-modal-panel cm-action-modal" data-modal="new-card" aria-hidden="true">',
    '    <div class="cm-modal-header">',
    '      <h2>New Card</h2>',
    '    </div>',
    '    <div class="cm-modal-body">',
           renderEditorForm(registry, null, "new"),
    '    </div>',
    '    <div class="cm-modal-actions">',
    '      <button class="cm-button cm-button-primary" type="button" data-modal-commit="new">Create</button>',
    '      <button class="cm-button cm-button-secondary" type="button" data-modal-close>Close</button>',
    '    </div>',
    '  </section>',
    '  <section class="cm-modal-panel cm-action-modal" data-modal="edit-card" aria-hidden="true">',
    '    <div class="cm-modal-header">',
    '      <h2>Edit Card</h2>',
    '    </div>',
    '    <div class="cm-modal-body" data-edit-card-body>',
    '      <div class="cm-modal-placeholder">Select a card to edit it.</div>',
    '    </div>',
    '    <div class="cm-modal-actions cm-modal-actions-split">',
    '      <div class="cm-modal-actions-left">',
    '        <button class="cm-button cm-button-secondary" type="button" data-delete-card disabled>Delete</button>',
    '      </div>',
    '      <div class="cm-modal-actions-right">',
    '        <button class="cm-button cm-button-primary" type="button" data-modal-commit="edit">Save</button>',
    '        <button class="cm-button cm-button-secondary" type="button" data-modal-close>Close</button>',
    '      </div>',
    '    </div>',
    '  </section>',
    '  <section class="cm-modal-panel cm-action-modal" data-modal="delete-card" aria-hidden="true">',
    '    <div class="cm-modal-header">',
    '      <h2>Delete Card</h2>',
    '    </div>',
    '    <div class="cm-modal-body">',
    '      <p class="cm-delete-copy">This will delete the selected card JSON file from your local repo. Git can still recover it if you reset or restore the file.</p>',
    '      <div class="cm-form-error" data-delete-error hidden></div>',
    '    </div>',
    '    <div class="cm-modal-actions">',
    '      <button class="cm-button cm-button-primary" type="button" data-confirm-delete-card>Delete</button>',
    '      <button class="cm-button cm-button-secondary" type="button" data-modal-close>Close</button>',
    '    </div>',
    '  </section>'
  ].join("");
}
