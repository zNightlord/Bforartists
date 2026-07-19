/* SPDX-FileCopyrightText: 2007 Blender Authors
 *
 * SPDX-License-Identifier: GPL-2.0-or-later */

/** \file
 * \ingroup nodes
 */

#pragma once

#include "BKE_node.hh"

namespace blender {

struct bNode;
struct bNodeTree;

bNodeSocket *node_add_socket_from_template(bNodeTree *ntree,
                                           bNode *node,
                                           bke::bNodeSocketTemplate *stemp,
                                           eNodeSocketInOut in_out);

void node_verify_sockets(Main *bmain, bNodeTree *ntree, bNode *node, bool do_id_user);

void node_socket_init_default_value_data(eNodeSocketDatatype datatype, int subtype, void **data);
void node_socket_copy_default_value_data(eNodeSocketDatatype datatype, void *to, const void *from);
/** Copy only the UI (subtype, soft range) of a default value struct, keeping #to's current value.
 * Used to refresh a socket's UI from a stored value struct without discarding the edited value. */
void node_socket_copy_default_value_ui_data(eNodeSocketDatatype datatype,
                                            void *to,
                                            const void *from);
/** Free a socket default value struct (e.g. #bNodeSocketValueFloat) for #datatype. Frees memory
 * only; ID user counts on ID-typed data are the caller's responsibility. */
void node_socket_free_default_value_data(eNodeSocketDatatype datatype, void *data);
/** Blend-write / read a socket default value struct for #datatype (the raw data; ID pointers are
 * resolved separately by lib linking). */
void node_socket_blend_write_default_value_data(BlendWriter *writer,
                                                eNodeSocketDatatype datatype,
                                                const void *data);
void node_socket_blend_read_default_value_data(BlendDataReader *reader,
                                               eNodeSocketDatatype datatype,
                                               void **data);
void node_socket_init_default_value(bNodeSocket *sock);
void node_socket_copy_default_value(bNodeSocket *to, const bNodeSocket *from);
void register_standard_node_socket_types();

namespace nodes {

/**
 * Change the sockets of the node so that it matches the declaration.
 *
 * \param bmain: Optional, necessary for updating animation data.
 */
void update_node_declaration_and_sockets(bNodeTree &ntree, bNode &node, Main *bmain = nullptr);
bool socket_type_supports_fields(eNodeSocketDatatype socket_type);
bool socket_type_supports_attributes(eNodeSocketDatatype socket_type);
bool socket_type_supports_grids(eNodeSocketDatatype socket_type);

}  // namespace nodes
}  // namespace blender
