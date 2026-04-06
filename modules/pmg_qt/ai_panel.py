"""
AI Chat Panel for PyMOL.

Provides a dockable QDockWidget with a chat interface that talks to an LLM
agent capable of running PyMOL commands and capturing the viewport.
"""

from pymol.Qt import QtWidgets, QtCore, QtGui

Qt = QtCore.Qt


# ---------------------------------------------------------------------------
# Worker thread for running the agent without blocking the UI
# ---------------------------------------------------------------------------

class _AgentWorker(QtCore.QThread):
    """Runs AIAgent.send_message() in a background thread."""

    finished = QtCore.Signal(str)       # emits final response text
    tool_called = QtCore.Signal(str)    # emits human-readable tool description
    reasoning = QtCore.Signal(str)      # emits model reasoning summary text
    error = QtCore.Signal(str)          # emits error message

    def __init__(self, agent, user_text: str, parent=None):
        super().__init__(parent)
        self._agent = agent
        self._user_text = user_text

    def run(self):
        try:
            response = self._agent.send_message(
                self._user_text,
                on_tool_call=lambda desc: self.tool_called.emit(desc),
                on_reasoning=lambda text: self.reasoning.emit(text),
            )
            self.finished.emit(response or "")
        except InterruptedError:
            self.finished.emit("")   # cancelled — emit empty, UI handles it
        except Exception as exc:  # noqa: BLE001
            self.error.emit(str(exc))


# ---------------------------------------------------------------------------
# Collapsible action log widget
# ---------------------------------------------------------------------------

class _CollapsibleActions(QtWidgets.QFrame):
    """
    A collapsible widget that shows a summary header (e.g. "3 actions taken")
    and expands to reveal individual tool-call entries. Collapsed by default.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._entries = []

        self.setStyleSheet(
            "QFrame { background-color: #2a2a2a; border-radius: 6px; border: 1px solid #3d6b3d; }"
        )

        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        # --- Header button (always visible, acts as toggle) ---
        self._toggle_btn = QtWidgets.QPushButton()
        self._toggle_btn.setFlat(True)
        self._toggle_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._toggle_btn.setStyleSheet(
            "QPushButton {"
            "  color: #6dbf6d;"
            "  font-size: 11px;"
            "  text-align: left;"
            "  padding: 5px 8px;"
            "  border: none;"
            "  background: transparent;"
            "}"
            "QPushButton:hover { color: #90d990; }"
        )
        self._toggle_btn.clicked.connect(self._toggle)
        outer.addWidget(self._toggle_btn)

        # --- Body (hidden by default) ---
        self._body = QtWidgets.QWidget()
        self._body.setVisible(False)
        body_layout = QtWidgets.QVBoxLayout(self._body)
        body_layout.setContentsMargins(8, 2, 8, 6)
        body_layout.setSpacing(2)
        self._body_layout = body_layout
        outer.addWidget(self._body)

        self._update_header()

    def add_action(self, description: str):
        """Append a new tool-call entry."""
        self._entries.append(description)

        label = QtWidgets.QLabel(f"• {description}")
        label.setWordWrap(True)
        label.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Ignored,
            QtWidgets.QSizePolicy.Policy.Preferred,
        )
        label.setStyleSheet(
            "QLabel { color: #8fbf8f; font-size: 11px; font-family: monospace; border: none; }"
        )
        self._body_layout.addWidget(label)
        self._update_header()

    def _update_header(self):
        n = len(self._entries)
        expanded = self._body.isVisible()
        arrow = "▼" if expanded else "▶"
        noun = "action" if n == 1 else "actions"
        self._toggle_btn.setText(f"{arrow}  {n} {noun} taken")

    def _toggle(self):
        visible = not self._body.isVisible()
        self._body.setVisible(visible)
        self._update_header()
        # Let the scroll area reflow
        self.updateGeometry()
        if parent := self.parent():
            parent.updateGeometry()


# ---------------------------------------------------------------------------
# Collapsible reasoning widget
# ---------------------------------------------------------------------------

class _CollapsibleReasoning(QtWidgets.QFrame):
    """
    A collapsible widget showing the model's intermediate reasoning summaries.
    Styled in a muted purple to distinguish from green action logs.
    Collapsed by default.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        self._chunks = []

        self.setStyleSheet(
            "QFrame { background-color: #251f2e; border-radius: 6px; border: 1px solid #5a3f7a; }"
        )

        outer = QtWidgets.QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self._toggle_btn = QtWidgets.QPushButton("▶  Thinking...")
        self._toggle_btn.setFlat(True)
        self._toggle_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._toggle_btn.setStyleSheet(
            "QPushButton {"
            "  color: #a07cc5;"
            "  font-size: 11px;"
            "  text-align: left;"
            "  padding: 5px 8px;"
            "  border: none;"
            "  background: transparent;"
            "}"
            "QPushButton:hover { color: #c4a0e8; }"
        )
        self._toggle_btn.clicked.connect(self._toggle)
        outer.addWidget(self._toggle_btn)

        self._body = QtWidgets.QWidget()
        self._body.setVisible(False)
        body_layout = QtWidgets.QVBoxLayout(self._body)
        body_layout.setContentsMargins(8, 2, 8, 6)
        body_layout.setSpacing(4)
        self._body_layout = body_layout
        outer.addWidget(self._body)

    def add_chunk(self, text: str):
        self._chunks.append(text)
        label = QtWidgets.QLabel(text)
        label.setWordWrap(True)
        label.setSizePolicy(
            QtWidgets.QSizePolicy.Policy.Ignored,
            QtWidgets.QSizePolicy.Policy.Preferred,
        )
        label.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse |
            Qt.TextInteractionFlag.TextSelectableByKeyboard
        )
        label.setStyleSheet(
            "QLabel { color: #b89fd4; font-size: 11px; font-style: italic; border: none; }"
        )
        self._body_layout.addWidget(label)

    def mark_done(self):
        """Update header once the agent has finished thinking."""
        arrow = "▼" if self._body.isVisible() else "▶"
        self._toggle_btn.setText(f"{arrow}  Thought ({len(self._chunks)} step(s))")

    def _toggle(self):
        visible = not self._body.isVisible()
        self._body.setVisible(visible)
        # Sync arrow — preserve "Thinking..." vs "Thought" label
        current = self._toggle_btn.text()
        if "▶" in current:
            self._toggle_btn.setText(current.replace("▶", "▼"))
        else:
            self._toggle_btn.setText(current.replace("▼", "▶"))
        self.updateGeometry()
        if parent := self.parent():
            parent.updateGeometry()


# ---------------------------------------------------------------------------
# Chat bubble widget
# ---------------------------------------------------------------------------

class _ChatBubble(QtWidgets.QFrame):
    """A single chat message rendered as a styled label."""

    def __init__(self, text: str, role: str, parent=None):
        """
        Args:
            text: Message content.
            role: "user" or "assistant" or "system"
        """
        super().__init__(parent)
        self.setMinimumWidth(0)

        label = QtWidgets.QLabel(text)
        label.setWordWrap(True)
        label.setMinimumWidth(0)
        label.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse |
            Qt.TextInteractionFlag.TextSelectableByKeyboard
        )

        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(8, 6, 8, 6)
        layout.addWidget(label)

        if role == "user":
            self.setStyleSheet(
                "QFrame { background-color: #2a5caa; border-radius: 8px; }"
                "QLabel { color: white; }"
            )
        elif role == "assistant":
            self.setStyleSheet(
                "QFrame { background-color: #3a3a3a; border-radius: 8px; }"
                "QLabel { color: #e0e0e0; }"
            )
        else:  # system / info
            self.setStyleSheet(
                "QFrame { background-color: #444; border-radius: 4px; }"
                "QLabel { color: #aaa; font-style: italic; font-size: 11px; }"
            )


# ---------------------------------------------------------------------------
# Main panel widget
# ---------------------------------------------------------------------------

class _AIPanelWidget(QtWidgets.QWidget):
    """Inner widget for the AI chat panel."""

    def __init__(self, cmd, parent=None):
        super().__init__(parent)
        self._cmd = cmd
        self._worker = None
        self._agent = None
        self._current_actions: "_CollapsibleActions | None" = None
        self._current_reasoning: "_CollapsibleReasoning | None" = None
        self._last_reasoning: "_CollapsibleReasoning | None" = None
        self._init_agent()
        self._build_ui()

    def _init_agent(self):
        try:
            from pmg_qt.ai_agent import AIAgent
            self._agent = AIAgent(self._cmd)
        except Exception as exc:  # noqa: BLE001
            self._agent = None
            self._init_error = str(exc)
        else:
            self._init_error = None

    def _build_ui(self):
        layout = QtWidgets.QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        # --- Toolbar ---
        toolbar = QtWidgets.QHBoxLayout()
        toolbar.setSpacing(6)

        reasoning_label = QtWidgets.QLabel("Reasoning:")
        reasoning_label.setStyleSheet("color: #aaa; font-size: 11px;")
        toolbar.addWidget(reasoning_label)

        self._reasoning_combo = QtWidgets.QComboBox()
        self._reasoning_combo.addItems(["high", "medium", "low", "none"])
        self._reasoning_combo.setCurrentText("high")
        self._reasoning_combo.setFixedWidth(80)
        self._reasoning_combo.setStyleSheet(
            "QComboBox { font-size: 11px; background-color: #2a2a2a; color: #ddd; border: 1px solid #555; border-radius: 3px; padding: 1px 4px; }"
            "QComboBox::drop-down { border: none; }"
            "QComboBox QAbstractItemView { background-color: #2a2a2a; color: #ddd; selection-background-color: #444; }"
        )
        self._reasoning_combo.currentTextChanged.connect(self._on_reasoning_mode_changed)
        toolbar.addWidget(self._reasoning_combo)
        toolbar.addStretch()

        layout.addLayout(toolbar)

        # --- Chat history area ---
        self._scroll_area = QtWidgets.QScrollArea()
        self._scroll_area.setWidgetResizable(True)
        self._scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

        self._history_widget = QtWidgets.QWidget()
        self._history_layout = QtWidgets.QVBoxLayout(self._history_widget)
        self._history_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        self._history_layout.setSpacing(6)
        self._history_layout.setContentsMargins(4, 4, 4, 4)

        self._scroll_area.setWidget(self._history_widget)
        layout.addWidget(self._scroll_area, stretch=1)

        # --- Status label ---
        self._status_label = QtWidgets.QLabel("")
        self._status_label.setStyleSheet("color: #aaa; font-style: italic; font-size: 11px;")
        layout.addWidget(self._status_label)

        # --- Input area ---
        input_layout = QtWidgets.QHBoxLayout()

        self._input = QtWidgets.QTextEdit()
        self._input.setFixedHeight(70)
        self._input.setPlaceholderText("Ask the AI anything about your structure...")
        self._input.installEventFilter(self)
        input_layout.addWidget(self._input, stretch=1)

        btn_layout = QtWidgets.QVBoxLayout()

        self._send_btn = QtWidgets.QPushButton("Send")
        self._send_btn.setFixedWidth(60)
        self._send_btn.clicked.connect(self._on_send)

        self._cancel_btn = QtWidgets.QPushButton("Cancel")
        self._cancel_btn.setFixedWidth(60)
        self._cancel_btn.setStyleSheet("QPushButton { color: #ff6b6b; }")
        self._cancel_btn.clicked.connect(self._on_cancel)
        self._cancel_btn.setVisible(False)

        self._clear_btn = QtWidgets.QPushButton("Clear")
        self._clear_btn.setFixedWidth(60)
        self._clear_btn.clicked.connect(self._on_clear)

        btn_layout.addWidget(self._send_btn)
        btn_layout.addWidget(self._cancel_btn)
        btn_layout.addWidget(self._clear_btn)
        btn_layout.addStretch()

        input_layout.addLayout(btn_layout)
        layout.addLayout(input_layout)

        # Show init error if agent failed to load
        if self._init_error:
            self._add_bubble(
                f"Error initializing AI agent: {self._init_error}\n\n"
                "Make sure 'openai' is installed (pip install openai) and "
                "OPENAI_API_KEY is set.",
                "system",
            )
            self._input.setEnabled(False)
            self._send_btn.setEnabled(False)

    # ------------------------------------------------------------------
    # Event filter: Ctrl+Return submits
    # ------------------------------------------------------------------

    def eventFilter(self, obj, event):
        if obj is self._input and event.type() == QtCore.QEvent.Type.KeyPress:
            if (event.key() in (Qt.Key.Key_Return, Qt.Key.Key_Enter) and
                    event.modifiers() & Qt.KeyboardModifier.ControlModifier):
                self._on_send()
                return True
        return super().eventFilter(obj, event)

    # ------------------------------------------------------------------
    # Slots
    # ------------------------------------------------------------------

    def _on_send(self):
        text = self._input.toPlainText().strip()
        if not text or self._worker is not None:
            return

        self._input.clear()
        self._add_bubble(text, "user")
        self._current_actions = None    # fresh action group for this turn
        self._current_reasoning = None  # fresh reasoning group for this turn
        self._last_reasoning = None
        self._set_running(True)
        self._status_label.setText("Thinking...")

        self._worker = _AgentWorker(self._agent, text, self)
        self._worker.finished.connect(self._on_agent_finished)
        self._worker.tool_called.connect(self._on_tool_called)
        self._worker.reasoning.connect(self._on_reasoning)
        self._worker.error.connect(self._on_agent_error)
        self._worker.start()

    def _on_reasoning_mode_changed(self, mode: str):
        if self._agent:
            self._agent.set_reasoning_effort(mode)

    def _on_cancel(self):
        if self._agent:
            self._agent.cancel()
        self._status_label.setText("Cancelling...")
        self._cancel_btn.setEnabled(False)

    def _on_reasoning(self, text: str):
        # Each new reasoning emission gets its own widget (interleaved with actions).
        # If actions were shown since the last reasoning, start a fresh reasoning widget.
        if self._current_reasoning is None:
            self._current_reasoning = _CollapsibleReasoning(self._history_widget)
            self._last_reasoning = self._current_reasoning
            self._history_layout.addWidget(self._current_reasoning)
        self._current_reasoning.add_chunk(text)
        # Reset actions so the next tool call starts a new action group after this reasoning
        self._current_actions = None
        QtCore.QTimer.singleShot(50, self._scroll_to_bottom)

    def _on_agent_finished(self, response: str):
        self._worker = None
        if self._last_reasoning is not None:
            self._last_reasoning.mark_done()
        self._current_actions = None
        self._current_reasoning = None
        self._last_reasoning = None
        self._status_label.setText("")
        self._set_running(False)
        if response:
            self._add_bubble(response, "assistant")
        else:
            self._add_bubble("Cancelled.", "system")
        self._input.setFocus()

    def _on_tool_called(self, description: str):
        self._status_label.setText(description)
        # Lazily create the collapsible action group. Reset reasoning so the next
        # reasoning chunk starts a fresh widget after this group of actions.
        if self._current_actions is None:
            self._current_actions = _CollapsibleActions(self._history_widget)
            self._history_layout.addWidget(self._current_actions)
        self._current_actions.add_action(description)
        # Reset reasoning so next reasoning chunk appears after these actions
        self._current_reasoning = None
        QtCore.QTimer.singleShot(50, self._scroll_to_bottom)

    def _on_agent_error(self, error_msg: str):
        self._worker = None
        self._current_actions = None
        self._current_reasoning = None
        self._last_reasoning = None
        self._status_label.setText("")
        self._set_running(False)
        self._add_bubble(f"Error: {error_msg}", "system")

    def _on_clear(self):
        self._current_actions = None
        self._current_reasoning = None
        self._last_reasoning = None
        if self._agent:
            self._agent.clear_history()
        while self._history_layout.count():
            item = self._history_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _add_bubble(self, text: str, role: str):
        bubble = _ChatBubble(text, role, self._history_widget)
        self._history_layout.addWidget(bubble)
        QtCore.QTimer.singleShot(50, self._scroll_to_bottom)

    def _scroll_to_bottom(self):
        vbar = self._scroll_area.verticalScrollBar()
        vbar.setValue(vbar.maximum())

    def _set_running(self, running: bool):
        """Toggle between idle state (Send visible) and running state (Cancel visible)."""
        self._input.setEnabled(not running)
        self._send_btn.setVisible(not running)
        self._cancel_btn.setVisible(running)
        self._cancel_btn.setEnabled(running)
        self._clear_btn.setEnabled(not running)
        self._reasoning_combo.setEnabled(not running)


# ---------------------------------------------------------------------------
# Public: QDockWidget wrapper (matches BuilderPanelDocked pattern)
# ---------------------------------------------------------------------------

class AIPanel(QtWidgets.QDockWidget):
    """Dockable AI chat panel for PyMOL."""

    def __init__(self, parent):
        super().__init__(parent)
        self.setWindowTitle("AI Assistant")
        self.setMinimumWidth(280)
        self.setMaximumWidth(500)

        inner = _AIPanelWidget(parent.cmd, self)
        self.setWidget(inner)
