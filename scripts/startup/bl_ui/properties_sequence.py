# SPDX-FileCopyrightText: 2018-2023 Blender Authors
#
# SPDX-License-Identifier: GPL-2.0-or-later

import bpy
from bpy.types import Menu, Panel, UIList

from bpy.app.translations import (
    contexts as i18n_contexts,
    pgettext_iface as iface_,
)

class SequenceButtonsPanel:
    bl_space_type = 'PROPERTIES'
    bl_region_type = 'WINDOW'
    bl_context = "sequence"
    # COMPAT_ENGINES must be defined in each subclass, external engines can add themselves here

    @classmethod
    def poll(cls, context):
        return (context.engine in cls.COMPAT_ENGINES)


class SEQUENCE_PT_format(SequenceButtonsPanel, Panel):
    bl_label = "Format"
    COMPAT_ENGINES = {
        'BLENDER_RENDER',
        'BLENDER_EEVEE_NEXT',
        'BLENDER_WORKBENCH',
    }

    _frame_rate_args_prev = None
    _preset_class = None

    @staticmethod
    def _draw_framerate_label(*args):
        # avoids re-creating text string each draw
        if SEQUENCE_PT_format._frame_rate_args_prev == args:
            return SEQUENCE_PT_format._frame_rate_ret

        fps, fps_base, preset_label = args

        if fps_base == 1.0:
            fps_rate = round(fps)
        else:
            fps_rate = round(fps / fps_base, 2)

        # TODO: Change the following to iterate over existing presets
        custom_framerate = (fps_rate not in {6, 8, 12, 23.98, 24, 25, 29.97, 30, 50, 59.94, 60, 120, 240})

        if custom_framerate is True:
            fps_label_text = iface_("Custom ({:.4g} fps)").format(fps_rate)
            show_framerate = True
        else:
            fps_label_text = iface_("{:.4g} fps").format(fps_rate)
            show_framerate = (preset_label == "Custom")

        SEQUENCE_PT_format._frame_rate_args_prev = args
        SEQUENCE_PT_format._frame_rate_ret = args = (fps_label_text, show_framerate)
        return args

    @staticmethod
    def draw_framerate(layout, rd):
        if SEQUENCE_PT_format._preset_class is None:
            SEQUENCE_PT_format._preset_class = bpy.types.RENDER_MT_framerate_presets

        args = rd.fps, rd.fps_base, SEQUENCE_PT_format._preset_class.bl_label
        fps_label_text, show_framerate = SEQUENCE_PT_format._draw_framerate_label(*args)

        if show_framerate:
            col = layout.column(align=True)
            col.prop(rd, "fps")
            col.prop(rd, "fps_base", text="Base")

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False  # No animation.

        rd = context.scene.render

        col = layout.column(align=True)
        col.prop(rd, "resolution_x", text="Resolution X")
        col.prop(rd, "resolution_y", text="Y")
        col.prop(rd, "resolution_percentage", text="%")

        col = layout.column(align=True)
        col.prop(rd, "pixel_aspect_x", text="Aspect X")
        col.prop(rd, "pixel_aspect_y", text="Y")

        col = layout.column(heading="Frame Rate")
        self.draw_framerate(col, rd)


class SEQUENCE_PT_frame_range(SequenceButtonsPanel, Panel):
    bl_label = "Frame Range"
    COMPAT_ENGINES = {
        'BLENDER_RENDER',
        'BLENDER_EEVEE_NEXT',
        'BLENDER_WORKBENCH',
    }

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False  # No animation.

        scene = context.scene

        col = layout.column(align=True)
        col.prop(scene, "frame_start", text="Frame Start")
        col.prop(scene, "frame_end", text="End")
        col.prop(scene, "frame_step", text="Step")


class SEQUENCE_PT_output(SequenceButtonsPanel, Panel):
    bl_label = "Output"
    COMPAT_ENGINES = {
        'BLENDER_RENDER',
        'BLENDER_EEVEE_NEXT',
        'BLENDER_WORKBENCH',
    }

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = False
        layout.use_property_decorate = False  # No animation.

        rd = context.scene.render
        image_settings = rd.image_settings

        layout.prop(rd, "filepath", text="")

        layout.use_property_split = True

        col = layout.column(heading="Saving")
        col.prop(rd, "use_file_extension")
        col.prop(rd, "use_render_cache")

        layout.template_image_settings(image_settings, color_management=False)

        if not rd.is_movie_format:
            col = layout.column(heading="Image Sequence")
            col.prop(rd, "use_overwrite")
            col.prop(rd, "use_placeholder")


class SEQUENCE_PT_output_views(SequenceButtonsPanel, Panel):
    bl_label = "Views"
    bl_parent_id = "SEQUENCE_PT_output"
    COMPAT_ENGINES = {
        'BLENDER_RENDER',
        'BLENDER_EEVEE_NEXT',
        'BLENDER_WORKBENCH',
    }

    @classmethod
    def poll(cls, context):
        rd = context.scene.render
        return rd.use_multiview

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = False
        layout.use_property_decorate = False  # No animation.

        rd = context.scene.render
        layout.template_image_views(rd.image_settings)


class SEQUENCE_PT_output_color_management(SequenceButtonsPanel, Panel):
    bl_label = "Color Management"
    bl_options = {'DEFAULT_CLOSED'}
    bl_parent_id = "SEQUENCE_PT_output"
    COMPAT_ENGINES = {
        'BLENDER_RENDER',
        'BLENDER_EEVEE_NEXT',
        'BLENDER_WORKBENCH',
    }

    def draw(self, context):
        image_settings = context.scene.render.image_settings

        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False  # No animation.

        flow = layout.grid_flow(row_major=True, columns=0, even_columns=False, even_rows=False, align=True)

        owner = image_settings

        col = flow.column()

        if image_settings.has_linear_colorspace:
            if hasattr(owner, "linear_colorspace_settings"):
                col.prop(owner.linear_colorspace_settings, "name", text="Color Space")
        else:
            col.prop(owner.display_settings, "display_device")
            col.separator()
            col.template_colormanaged_view_settings(owner, "view_settings")


class SEQUENCE_PT_encoding(SequenceButtonsPanel, Panel):
    bl_label = "Encoding"
    bl_parent_id = "SEQUENCE_PT_output"
    bl_options = {'DEFAULT_CLOSED'}
    COMPAT_ENGINES = {
        'BLENDER_RENDER',
        'BLENDER_EEVEE_NEXT',
        'BLENDER_WORKBENCH',
    }

    @classmethod
    def poll(cls, context):
        rd = context.scene.render
        return rd.image_settings.file_format in {'FFMPEG', 'XVID', 'H264', 'THEORA'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False

        rd = context.scene.render
        ffmpeg = rd.ffmpeg

        layout.prop(rd.ffmpeg, "format")
        layout.prop(ffmpeg, "use_autosplit")


class SEQUENCE_PT_encoding_video(SequenceButtonsPanel, Panel):
    bl_label = "Video"
    bl_parent_id = "SEQUENCE_PT_encoding"
    COMPAT_ENGINES = {
        'BLENDER_RENDER',
        'BLENDER_EEVEE_NEXT',
        'BLENDER_WORKBENCH',
    }

    @classmethod
    def poll(cls, context):
        rd = context.scene.render
        return rd.image_settings.file_format in {'FFMPEG', 'XVID', 'H264', 'THEORA'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False

        self.draw_vcodec(context)

    def draw_vcodec(self, context):
        """Video codec options."""
        layout = self.layout
        ffmpeg = context.scene.render.ffmpeg

        needs_codec = ffmpeg.format in {
            'AVI',
            'QUICKTIME',
            'MKV',
            'OGG',
            'MPEG4',
            'WEBM',
        }
        if needs_codec:
            layout.prop(ffmpeg, "codec")

        if needs_codec and ffmpeg.codec == 'NONE':
            return

        # Color depth. List of codecs needs to be in sync with
        # `IMB_ffmpeg_valid_bit_depths` in source code.
        use_bpp = needs_codec and ffmpeg.codec in {'H264', 'H265', 'AV1', 'PRORES'}
        if use_bpp:
            image_settings = context.scene.render.image_settings
            layout.prop(image_settings, "color_depth", expand=True)

        if ffmpeg.codec == 'DNXHD':
            layout.prop(ffmpeg, "use_lossless_output")

        if ffmpeg.codec == 'PRORES':
            layout.prop(ffmpeg, "ffmpeg_prores_profile")

        # Output quality
        use_crf = needs_codec and ffmpeg.codec in {
            'H264',
            'H265',
            'MPEG4',
            'WEBM',
            'AV1',
        }
        if use_crf:
            layout.prop(ffmpeg, "constant_rate_factor")

        use_encoding_speed = needs_codec and ffmpeg.codec not in {'DNXHD', 'FFV1', 'HUFFYUV', 'PNG', 'PRORES', 'QTRLE'}
        use_bitrate = needs_codec and ffmpeg.codec not in {'FFV1', 'HUFFYUV', 'PNG', 'PRORES', 'QTRLE'}
        use_min_max_bitrate = ffmpeg.codec not in {'DNXHD'}
        use_gop = needs_codec and ffmpeg.codec not in {'DNXHD', 'HUFFYUV', 'PNG', 'PRORES'}
        use_b_frames = needs_codec and use_gop and ffmpeg.codec not in {'FFV1', 'QTRLE'}

        # Encoding speed
        if use_encoding_speed:
            layout.prop(ffmpeg, "ffmpeg_preset")
        # I-frames
        if use_gop:
            layout.prop(ffmpeg, "gopsize")
        # B-Frames
        if use_b_frames:
            row = layout.row(align=True, heading="Max B-frames")
            row.prop(ffmpeg, "use_max_b_frames", text="")
            sub = row.row(align=True)
            sub.active = ffmpeg.use_max_b_frames
            sub.prop(ffmpeg, "max_b_frames", text="")

        if (not use_crf or ffmpeg.constant_rate_factor == 'NONE') and use_bitrate:
            col = layout.column()

            sub = col.column(align=True)
            sub.prop(ffmpeg, "video_bitrate")
            if use_min_max_bitrate:
                sub.prop(ffmpeg, "minrate", text="Minimum")
                sub.prop(ffmpeg, "maxrate", text="Maximum")

                col.prop(ffmpeg, "buffersize", text="Buffer")

                col.separator()

                col.prop(ffmpeg, "muxrate", text="Mux Rate")
                col.prop(ffmpeg, "packetsize", text="Mux Packet Size")


class SEQUENCE_PT_encoding_audio(SequenceButtonsPanel, Panel):
    bl_label = "Audio"
    bl_parent_id = "SEQUENCE_PT_encoding"
    COMPAT_ENGINES = {
        'BLENDER_RENDER',
        'BLENDER_EEVEE_NEXT',
        'BLENDER_WORKBENCH',
    }

    @classmethod
    def poll(cls, context):
        rd = context.scene.render
        return rd.image_settings.file_format in {'FFMPEG', 'XVID', 'H264', 'THEORA'}

    def draw(self, context):
        layout = self.layout
        layout.use_property_split = True
        layout.use_property_decorate = False

        rd = context.scene.render
        ffmpeg = rd.ffmpeg

        if ffmpeg.format != 'MP3':
            layout.prop(ffmpeg, "audio_codec", text="Audio Codec")

        if ffmpeg.audio_codec != 'NONE':
            layout.prop(ffmpeg, "audio_channels")
            layout.prop(ffmpeg, "audio_mixrate", text="Sample Rate")
            layout.prop(ffmpeg, "audio_bitrate")
            layout.prop(ffmpeg, "audio_volume", slider=True)

classes = (
    SEQUENCE_PT_format,
    SEQUENCE_PT_frame_range,
    SEQUENCE_PT_output,
    SEQUENCE_PT_output_views,
    SEQUENCE_PT_output_color_management,
    SEQUENCE_PT_encoding,
    SEQUENCE_PT_encoding_video,
    SEQUENCE_PT_encoding_audio,
)

if __name__ == "__main__":  # only for live edit.
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
