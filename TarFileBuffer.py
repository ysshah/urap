# ===================== A utility to handle tar file buffer instead of making many small files =======
class TarFileBuffer:
    class ArchiveFile:
        def __init__(self, name, parent_tarfile):
            import io
            self.buffer = io.BytesIO()
            self.buffer.name = name
            self.parent_tarfile = parent_tarfile

        def close(self):
            # Flush the buffer and register the content to the tar file.
            import tarfile, time
            self.buffer.flush()
            self.buffer.seek(0)
            buf_values = self.buffer.getvalue()
            tarinfo = tarfile.TarInfo(name=self.buffer.name)
            tarinfo.size = len(buf_values)
            tarinfo.mtime = time.time()
            self.parent_tarfile.addfile(tarinfo=tarinfo, fileobj=self.buffer)
            self.buffer.close()
            self.buffer = None

        def __del__(self):
            if self.buffer != None:
                raise RuntimeError('TarFileBuffer.ArchiveFile: You probably didn\'t call close().')

    def __init__(self, name):
        import tarfile
        self.filebuf = tarfile.TarFile(name,'w')

    def newfile(self, name):
        return self.ArchiveFile(name, self.filebuf)

    def close(self):
        self.filebuf.close()
        self.filebuf = None

    def __del__(self):
        if self.filebuf!=None:
            self.close()
# =====================================================
